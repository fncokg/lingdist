#include <Rcpp.h>
#include <RcppThread.h>

#include "helpers.hpp"
#include "indexed_costs.hpp"
#include "prepare_input.hpp"
#include "cost_io.hpp"
#include "dp.hpp"
#include "expected_counts.hpp"
#include "dist.hpp"
#include "parallel.hpp"
#include <algorithm>
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
            if (a < b)
                std::swap(a, b);
            return a * (a + 1) / 2 + b;
        }

        void add(lingdist::TokenId a, lingdist::TokenId b, double weight)
        {
            if (a == b)
                return;
            pair_counts[index(a, b)] += weight;
            sym_counts[a] += weight;
            sym_counts[b] += weight;
        }

        void merge(const RawAlignmentCounts &other)
        {
            for (std::size_t i = 0; i < pair_counts.size(); ++i)
                pair_counts[i] += other.pair_counts[i];
            for (std::size_t i = 0; i < sym_counts.size(); ++i)
                sym_counts[i] += other.sym_counts[i];
        }
    };

    struct TrainingBlock
    {
        lingdist::ExpectedCountWorkspace workspace;
        RawAlignmentCounts counts;

        explicit TrainingBlock(std::size_t symbols) : counts(symbols) {}
    };

    lingdist::IndexedCosts update_costs(const RawAlignmentCounts &counts, const lingdist::IndexedCosts &previous,
                                        double damping_alpha, double log_smoothing_lambda, double min_edit_cost)
    {
        if (damping_alpha == 0.0 || counts.nsyms < 2)
            return previous;
        const auto &symbols = previous.symbol_table();
        // Add lambda to each unordered nonmatching pair, including character-gap pairs but excluding every diagonal.

        // Deal with marginal gains once: Each symbol is an endpoint in V-1 such pairs, so its marginal gains (V-1)*lambda, i.e. log(V-1)+log(lambda) in log space.
        const double log_marginal_prior = std::log(static_cast<double>(counts.nsyms - 1)) + log_smoothing_lambda;
        std::vector<double> log_marginals(counts.nsyms);
        for (std::size_t i = 0; i < counts.nsyms; ++i)
            log_marginals[i] = lingdist::logaddexp(std::log(counts.sym_counts[i]), log_marginal_prior);

        std::vector<double> scores(counts.pair_counts.size(), 0.0);
        double min_raw = std::numeric_limits<double>::infinity();
        double max_raw = -std::numeric_limits<double>::infinity();
        for (std::size_t r = 0; r < counts.nsyms; ++r)
            for (std::size_t c = 0; c < r; ++c)
            {
                const auto i = counts.index(r, c);
                // Total counts serve as denominators for both join and marginal probabilities. They can be omitted because they cancel out in the latter max-min normalization since we are in log space.
                // For the same reason, the base of the logarithm is irrelevant.

                // $\mathrm{PMI}=-\log\frac{p(x,y)}{p(x)p(y)}$

                // here we directly compute $\log c(x)+\log c(y)-\log c(x,y)$, of course, smoothed version.

                const double log_joint = lingdist::logaddexp(std::log(counts.pair_counts[i]), log_smoothing_lambda);
                const double raw = log_marginals[r] + log_marginals[c] - log_joint;
                min_raw = std::min(min_raw, raw);
                max_raw = std::max(max_raw, raw);
                scores[i] = raw;
            }

        // The candidate array is mutable only while being constructed. Moving it into a new snapshot leaves the previous epoch and its caches valid.
        std::vector<double> values(counts.nsyms * counts.nsyms, 0.0);
        for (std::size_t r = 0; r < counts.nsyms; ++r)
            for (std::size_t c = 0; c < r; ++c)
            {
                // Shift the nonmatching target range to [min_edit_cost, 1].
                // The zero-filled diagonal remains 0
                const double target = lingdist::double_equal(min_raw, max_raw)
                                          ? 1.0
                                          : min_edit_cost + (1.0 - min_edit_cost) *
                                                                ((scores[counts.index(r, c)] - min_raw) / (max_raw - min_raw));
                // Although the damping formula preserves this range mathematically;
                // we need to clamp it to [min_edit_cost, 1] to avoid possible rounding errors, so as to retain the strictly positive lower bound.
                const double value = std::clamp(damping_alpha == 1.0 ? target
                                                                     : (1.0 - damping_alpha) * previous.get(r, c) + damping_alpha * target,
                                                min_edit_cost, 1.0);
                values[c * counts.nsyms + r] = value;
                values[r * counts.nsyms + c] = value;
            }
        return lingdist::IndexedCosts(symbols, std::move(values));
    }

}

namespace lingdist
{
    List pmi_df(const DataFrame &data, const String &delim, bool detailed, const String &normalize_method,
                bool squareform, int n_threads, int max_epochs, double tol, double damping_alpha, double smoothing_lambda,
                double min_edit_cost,
                const String &form_strategy, const String &form_delim, std::vector<double> form_weights, bool quiet)
    {
        if (max_epochs < 1)
            stop("max_epochs must be at least 1.");
        if (!std::isfinite(tol) || tol <= 0.0)
            stop("tol must be finite and positive.");
        if (!std::isfinite(damping_alpha) || damping_alpha < 0.0 || damping_alpha > 1.0)
            stop("damping_alpha must be finite and between 0 and 1 inclusive.");
        if (!std::isfinite(smoothing_lambda) || smoothing_lambda <= 0.0)
            stop("smoothing_lambda must be finite and strictly positive.");
        if (!std::isfinite(min_edit_cost) || min_edit_cost <= 0.0 || min_edit_cost > 1.0)
            stop("min_edit_cost must be finite and in (0, 1].");
        // Lambda is fixed across epochs. Keep its linear and log meanings distinct.
        const double log_smoothing_lambda = std::log(smoothing_lambda);
        if (n_threads < 1)
            stop("n_threads must be at least 1.");
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

        // Here in PMI, different with the standard edit distance, we divide all tasks into several blocks, instead of parallel for all tasks.
        // This is because we may need to accumulate the same `sym_counts` after every DP alignment, which is not thread-safe.
        // So we divide all tasks into several blocks, and each block will accumulate its own `sym_counts`, and then we merge them after all blocks are done.
        const auto available = std::max(1u, std::thread::hardware_concurrency());
        const auto block_count = std::min({plan.tasks.size(), static_cast<std::size_t>(n_threads),
                                           static_cast<std::size_t>(available)});
        std::vector<TrainingBlock> blocks;
        blocks.reserve(block_count);
        for (std::size_t i = 0; i < block_count; ++i)
            blocks.emplace_back(symbols.size());

        if (!quiet)
            Rprintf("Starting PMI iterations with %llu weighted word-pair tasks and %d symbols.\n",
                    static_cast<unsigned long long>(plan.tasks.size()), static_cast<int>(symbols.size()));
        double sum_diff = 0.0, mean_diff = 0.0;
        int n_epochs = 0;
        bool converged = false;
        for (int epoch = 1; epoch <= max_epochs; ++epoch)
        {
            if (!quiet)
                Rprintf("Epoch %d/%d\n", epoch, max_epochs);
            std::unique_ptr<SafeProgressBar> bar;
            if (!quiet)
                bar = std::make_unique<SafeProgressBar>(plan.tasks.size(), 1);
            parallel_for_checked(block_count, [&](std::size_t block_id)
                                 {
                auto &block = blocks[block_id];
                block.counts.clear();
                // Each numbered block owns its counts and scratch arrays,
                // regardless of which pool thread happens to execute it.
                const auto width = plan.tasks.size() / block_count;
                const auto remainder = plan.tasks.size() % block_count;
                const auto begin = block_id * width + std::min(block_id, remainder);
                const auto end = begin + width + (block_id < remainder ? 1 : 0);
                for (auto i = begin; i < end; ++i)
                {
                    // One parallel iteration owns a whole block, so also poll
                    // between its word-pair tasks for responsive cancellation.
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
                } }, n_threads);

            counts.clear();
            for (const auto &block : blocks)
            {
                counts.merge(block.counts); // Stable block order within this call.
            }
            const auto previous = cost; // Shares immutable values, no V-by-V copy.
            cost = update_costs(counts, previous, damping_alpha, log_smoothing_lambda, min_edit_cost);
            sum_diff = 0.0;
            for (std::size_t target = 0; target < symbols.size(); ++target)
                for (std::size_t source = 0; source < symbols.size(); ++source)
                    sum_diff += std::fabs(cost.get(source, target) - previous.get(source, target));
            const auto off_diagonal = symbols.size() * (symbols.size() - 1);
            mean_diff = off_diagonal ? sum_diff / off_diagonal : 0.0;
            n_epochs = epoch;
            if (!quiet)
                Rprintf("Cost snapshot updated, mean absolute difference: %.6f\n", mean_diff);
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
