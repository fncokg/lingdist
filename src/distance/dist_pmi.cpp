#include <Rcpp.h>
#include <RcppThread.h>

#include "helpers.hpp"
#include "indexed_costs.hpp"
#include "prepare_input.hpp"
#include "cost_io.hpp"
#include "dp.hpp"
#include "expected_counts.hpp"
#include "dist.hpp"
#include <algorithm>
#include <exception>
#include <thread>

using namespace Rcpp;

namespace
{

    struct RawAlignmentCounts
    {
        std::size_t nsyms;
        std::vector<double> pair_counts;
        std::vector<double> sym_counts;

        explicit RawAlignmentCounts(std::size_t size)
            : nsyms(size), pair_counts(size * (size + 1) / 2, 0.0), sym_counts(size, 0.0) {}

        void clear()
        {
            // Block-local buffers contain observations only. Adding priors here
            // would make their strength depend on the number of worker blocks.
            std::fill(pair_counts.begin(), pair_counts.end(), 0.0);
            std::fill(sym_counts.begin(), sym_counts.end(), 0.0);
        }

        std::size_t index(std::size_t a, std::size_t b) const
        {
            if (a < b) std::swap(a, b);
            return a * (a + 1) / 2 + b;
        }

        void add(lingdist::TokenId a, lingdist::TokenId b, double weight)
        {
            if (a == b) return;
            pair_counts[index(a, b)] += weight;
            sym_counts[a] += weight;
            sym_counts[b] += weight;
        }

        void merge(const RawAlignmentCounts &other)
        {
            for (std::size_t i = 0; i < pair_counts.size(); ++i) pair_counts[i] += other.pair_counts[i];
            for (std::size_t i = 0; i < sym_counts.size(); ++i) sym_counts[i] += other.sym_counts[i];
        }
    };

    struct TrainingBlock
    {
        lingdist::ExpectedCountWorkspace workspace;
        RawAlignmentCounts counts;
        std::exception_ptr error;

        explicit TrainingBlock(std::size_t symbols) : counts(symbols) {}
    };

    lingdist::IndexedCosts update_costs(const RawAlignmentCounts &counts, const lingdist::IndexedCosts &previous,
                                       double damping_alpha)
    {
        if (damping_alpha == 0.0) return previous;
        const auto &symbols = previous.symbol_table();
        // Apply the existing independent +1 pair/marginal priors once, after
        // merging all observations. Consistent smoothing is a separate change.
        // Keep raw counts intact rather than overwriting them with PMI scores.
        std::vector<double> scores(counts.pair_counts.size(), 0.0);
        double min_raw = std::numeric_limits<double>::infinity();
        double max_raw = -std::numeric_limits<double>::infinity();
        for (std::size_t r = 0; r < counts.nsyms; ++r)
            for (std::size_t c = 0; c < r; ++c)
            {
                const auto i = counts.index(r, c);
                // Constants common to all negative PMI entries disappear in
                // the subsequent min-max mapping; use the same formula as before.
                const double raw = -std::log((counts.pair_counts[i] + 1.0) /
                    ((counts.sym_counts[r] + 1.0) * (counts.sym_counts[c] + 1.0)));
                min_raw = std::min(min_raw, raw);
                max_raw = std::max(max_raw, raw);
                scores[i] = raw;
            }

        // The candidate array is mutable only while being constructed. Moving
        // it into a new snapshot leaves the previous epoch and its caches valid.
        std::vector<double> values(counts.nsyms * counts.nsyms, 0.0);
        for (std::size_t r = 0; r < counts.nsyms; ++r)
            for (std::size_t c = 0; c < r; ++c)
            {
                const double target = lingdist::double_equal(min_raw, max_raw)
                    ? 1.0 : (scores[counts.index(r, c)] - min_raw) / (max_raw - min_raw);
                // Damp the mapped target, so convergence and prediction use
                // W_next = (1-alpha)*W_previous + alpha*W_target. Keep the
                // alpha=1 path identical to the original full update.
                const double value = damping_alpha == 1.0 ? target
                    : (1.0 - damping_alpha) * previous.get(r, c) + damping_alpha * target;
                values[c * counts.nsyms + r] = value;
                values[r * counts.nsyms + c] = value;
            }
        return lingdist::IndexedCosts(symbols, std::move(values));
    }


}

namespace lingdist
{
    List pmi_df(const DataFrame &data, const String &delim, bool detailed, const String &normalize_method,
                bool squareform, int n_threads, int max_epochs, double tol, double damping_alpha,
                const String &form_strategy, const String &form_delim, std::vector<double> form_weights, bool quiet)
    {
        if (max_epochs < 1) stop("max_epochs must be at least 1.");
        if (!std::isfinite(tol) || tol <= 0.0) stop("tol must be finite and positive.");
        if (!std::isfinite(damping_alpha) || damping_alpha < 0.0 || damping_alpha > 1.0)
            stop("damping_alpha must be finite and between 0 and 1 inclusive.");
        if (n_threads < 1) stop("n_threads must be at least 1.");
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
        RawAlignmentCounts counts(symbols.size());
        const auto available = std::max(1u, std::thread::hardware_concurrency());
        const auto block_count = std::min({plan.tasks.size(), static_cast<std::size_t>(n_threads),
                                           static_cast<std::size_t>(available)});
        std::vector<TrainingBlock> blocks;
        blocks.reserve(block_count);
        for (std::size_t i = 0; i < block_count; ++i) blocks.emplace_back(symbols.size());
        if (!quiet)
            Rprintf("Starting PMI iterations with %llu weighted word-pair tasks and %d symbols.\n",
                    static_cast<unsigned long long>(plan.tasks.size()), static_cast<int>(symbols.size()));
        double sum_diff = 0.0, mean_diff = 0.0;
        int n_epochs = 0;
        bool converged = false;
        for (int epoch = 1; epoch <= max_epochs; ++epoch)
        {
            if (!quiet) Rprintf("Epoch %d/%d\n", epoch, max_epochs);
            std::unique_ptr<SafeProgressBar> bar;
            if (!quiet) bar = std::make_unique<SafeProgressBar>(plan.tasks.size(), 1);
            RcppThread::parallelFor(0, static_cast<int>(block_count), [&](int block_id)
            {
                auto &block = blocks[block_id];
                block.counts.clear();
                block.error = nullptr;
                // Each numbered block owns its counts and scratch arrays,
                // regardless of which pool thread happens to execute it.
                const auto width = plan.tasks.size() / block_count;
                const auto remainder = plan.tasks.size() % block_count;
                const auto id = static_cast<std::size_t>(block_id);
                const auto begin = id * width + std::min(id, remainder);
                const auto end = begin + width + (id < remainder ? 1 : 0);
                try
                {
                    for (auto i = begin; i < end; ++i)
                    {
                        if (RcppThread::isInterrupted()) return;
                        const auto &task = plan.tasks[i];
                        accumulate_expected_alignment_counts(input.word(task.pair.source), input.word(task.pair.target),
                            cost, block.workspace, [&](TokenId a, TokenId b, double probability)
                            {
                                // Form weights are already normalized per cell.
                                // Path averaging is entirely inside probability.
                                block.counts.add(a, b, task.weight * probability);
                            });
                        if (bar) (*bar)++;
                    }
                }
                catch (...)
                {
                    // C++ kernel failures must reach Rcpp on the calling thread,
                    // rather than escape a pool worker and terminate R.
                    block.error = std::current_exception();
                }
            }, n_threads);

            counts.clear();
            for (const auto &block : blocks)
            {
                if (block.error) std::rethrow_exception(block.error);
                counts.merge(block.counts); // Stable block order within this call.
            }
            const auto previous = cost; // Shares immutable values, no V-by-V copy.
            cost = update_costs(counts, previous, damping_alpha);
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
