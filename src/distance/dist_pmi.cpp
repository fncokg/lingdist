#include <Rcpp.h>
#include <RcppThread.h>

#include "helpers.hpp"
#include "indexed_costs.hpp"
#include "prepare_input.hpp"
#include "cost_io.hpp"
#include "dp.hpp"
#include "alignment.hpp"
#include "dist.hpp"

using namespace Rcpp;

namespace
{

    struct AlignmentCount
    {
        std::size_t nsyms;
        std::vector<double> pair_counts;
        std::vector<double> sym_counts;

        explicit AlignmentCount(std::size_t size)
            : nsyms(size), pair_counts(size * (size + 1) / 2, 1.0), sym_counts(size, 1.0) {}

        void clear()
        {
            // Retain the current independent priors. Consistent smoothing is a
            // separate statistical change, not part of replacing cost storage.
            std::fill(pair_counts.begin(), pair_counts.end(), 1.0);
            std::fill(sym_counts.begin(), sym_counts.end(), 1.0);
        }

        std::size_t index(std::size_t a, std::size_t b) const
        {
            if (a < b) std::swap(a, b);
            return a * (a + 1) / 2 + b;
        }

        void add(lingdist::TokenId a, lingdist::TokenId b, double weight)
        {
            pair_counts[index(a, b)] += weight;
            sym_counts[a] += weight;
            sym_counts[b] += weight;
        }
    };

    lingdist::IndexedCosts update_costs(const std::vector<lingdist::AlignmentResult> &results, const lingdist::TrainingPlan &plan,
                                       const lingdist::SymbolTable &symbols, AlignmentCount &counts)
    {
        counts.clear();
        for (std::size_t i = 0; i < results.size(); ++i)
        {
            const auto &word = results[i];
            if (word.alignments.empty()) continue;
            // Two distinct averages: the task weight sums normalized form-pair
            // weights across cells; 1/K averages this word pair's retained paths.
            const double weight = plan.tasks[i].weight / static_cast<double>(word.alignments.size());
            for (const auto &path : word.alignments)
                for (const auto &pair : path)
                    if (pair.first != pair.second) counts.add(pair.first, pair.second, weight);
        }

        double min_raw = std::numeric_limits<double>::infinity();
        double max_raw = -std::numeric_limits<double>::infinity();
        for (std::size_t r = 0; r < counts.nsyms; ++r)
            for (std::size_t c = 0; c < r; ++c)
            {
                const auto i = counts.index(r, c);
                // Constants common to all negative PMI entries disappear in
                // the subsequent min-max mapping; use the same formula as before.
                const double raw = -std::log(counts.pair_counts[i] / (counts.sym_counts[r] * counts.sym_counts[c]));
                min_raw = std::min(min_raw, raw);
                max_raw = std::max(max_raw, raw);
                counts.pair_counts[i] = raw;
            }

        // The candidate array is mutable only while being constructed. Moving
        // it into a new snapshot leaves the previous epoch and its caches valid.
        std::vector<double> values(counts.nsyms * counts.nsyms, 0.0);
        for (std::size_t r = 0; r < counts.nsyms; ++r)
            for (std::size_t c = 0; c < r; ++c)
            {
                const double value = lingdist::double_equal(min_raw, max_raw)
                    ? 1.0 : (counts.pair_counts[counts.index(r, c)] - min_raw) / (max_raw - min_raw);
                values[c * counts.nsyms + r] = value;
                values[r * counts.nsyms + c] = value;
            }
        return lingdist::IndexedCosts(symbols, std::move(values));
    }


}

namespace lingdist
{
    List pmi_df(const DataFrame &data, const String &delim, bool detailed, const String &normalize_method,
                bool squareform, int n_threads, int max_epochs, double tol, int alignment_max_paths,
                const String &form_strategy, const String &form_delim, std::vector<double> form_weights, bool quiet)
    {
        if (max_epochs < 1) stop("max_epochs must be at least 1.");
        if (!std::isfinite(tol) || tol <= 0.0) stop("tol must be finite and positive.");
        if (alignment_max_paths < 1) stop("alignment_max_paths must be at least 1.");
        if (normalize_method != "longest" && normalize_method != "none")
            stop("Unsupported normalize_method. Supported values are 'longest' and 'none'.");

        // Encode once, sharing the same vocabulary through every cost snapshot
        // and final prediction. Costs do not retain the word table or cell map.
        const auto strategy = parse_form_strategy(static_cast<std::string>(form_strategy));
        if (strategy == FormStrategy::Minimum)
            stop("PMI training supports weighting, mean_all, first or off; min_all depends on the current costs.");
        const auto input = prepare_input(data, delim, form_delim, strategy != FormStrategy::Off);
        const auto policy = make_form_policy(strategy, input.max_forms(), std::move(form_weights));
        const auto plan = build_training_plan(input, policy);
        const auto &symbols = input.symbol_table();
        auto cost = IndexedCosts::uniform(symbols);
        AlignmentCount counts(symbols.size());
        if (plan.tasks.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
            stop("Too many unique word-pair tasks for parallel execution.");
        const int n_tasks = static_cast<int>(plan.tasks.size());
        if (!quiet)
            Rprintf("Starting PMI iterations with %d weighted word-pair tasks and %d symbols.\n", n_tasks, static_cast<int>(symbols.size()));
        std::vector<AlignmentResult> results(plan.tasks.size());
        double sum_diff = 0.0, mean_diff = 0.0;
        int n_epochs = 0;
        bool converged = false;
        for (int epoch = 1; epoch <= max_epochs; ++epoch)
        {
            if (!quiet) Rprintf("Epoch %d/%d\n", epoch, max_epochs);
            std::unique_ptr<SafeProgressBar> bar;
            if (!quiet) bar = std::make_unique<SafeProgressBar>(plan.tasks.size(), 1);
            RcppThread::parallelFor(0, n_tasks, [&](int i)
            {
                const auto pair = plan.tasks[i].pair;
                results[i] = get_alignment_result(input.word(pair.source), input.word(pair.target), cost, alignment_max_paths, "none");
                if (bar) (*bar)++;
            }, n_threads);

            const auto previous = cost; // Shares immutable values, no V-by-V copy.
            cost = update_costs(results, plan, symbols, counts);
            sum_diff = 0.0;
            for (std::size_t target = 0; target < symbols.size(); ++target)
                for (std::size_t source = 0; source < symbols.size(); ++source)
                    sum_diff += std::fabs(cost.get(source, target) - previous.get(source, target));
            const auto off_diagonal = symbols.size() * (symbols.size() - 1);
            mean_diff = off_diagonal ? sum_diff / off_diagonal : 0.0;
            n_epochs = epoch;
            if (!quiet) Rprintf("Cost snapshot updated, mean absolute difference: %.6f\n", mean_diff);
            if (mean_diff < tol)
            {
                converged = true;
                break;
            }
        }
        if (!converged)
            warning("PMI distance computation did not converge within the maximum number of epochs. Try increasing max_epochs or tol.");

        // Evaluate with the returned snapshot, not the previous epoch's costs.
        const auto result = edit_dist_prepared(input, cost, plan.policy, detailed, normalize_method, squareform, true, n_threads, quiet);
        return List::create(Named("result") = result, Named("cost") = write_costs(cost),
                            Named("sum_diff") = sum_diff, Named("mean_diff") = mean_diff,
                            Named("converged") = converged, Named("n_epochs") = n_epochs);
    }
}
