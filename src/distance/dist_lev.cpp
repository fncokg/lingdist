#include <Rcpp.h>
#include <RcppThread.h>

#include "dist.hpp"
#include "dp.hpp"
#include "dispatcher.hpp"
#include "cost_io.hpp"
#include "prepare_input.hpp"
#include "parallel.hpp"
#include <unordered_map>

namespace
{
    void validate_normalization(const String &method)
    {
        if (method != "longest" && method != "none")
            stop("Unsupported normalize_method. Supported values are 'longest' and 'none'.");
    }

    void warn_missing(const lingdist::SymbolTable &symbols, const DataFrame &matrix)
    {
        const auto missing = lingdist::missing_cost_symbols(symbols, matrix);
        if (missing.empty())
            return;
        std::string message;
        for (const auto &name : missing)
            message += name + ", ";
        warning("The cost matrix is missing symbols: %s. Undefined operations use 0 for matching and 1 otherwise.", message.c_str());
    }

    struct DistancePlan
    {
        std::unordered_map<lingdist::WordPairKey, std::size_t, lingdist::WordPairHash> task_ids;
        std::vector<lingdist::WordPairKey> tasks;
    };

    DistancePlan build_distance_plan(const lingdist::PreparedInput &input,
                                     const std::vector<lingdist::RowPair> &pairs,
                                     const lingdist::FormPolicy &policy)
    {
        DistancePlan plan;
        for (std::size_t i = 0; i < pairs.size(); ++i)
        {
            if (i % 1024 == 0)
                Rcpp::checkUserInterrupt();
            for (std::size_t col = 0; col < input.ncol(); ++col)
                lingdist::for_each_form_pair(input.cell_forms(pairs[i].source, col),
                                             input.cell_forms(pairs[i].target, col), policy,
                                             [&](lingdist::WordId a, lingdist::WordId b, double)
                                             {
                    // Enumerate exactly the pairs used by prediction, including min_all.
                    // Keep direction even for symmetric output: that flag does not ensure the symmetry of the supplied edit costs.
                    const auto key = lingdist::make_word_pair_key(a, b);
                    const auto inserted = plan.task_ids.try_emplace(key, plan.tasks.size()).second;
                    if (inserted) plan.tasks.push_back(key); });
        }
        return plan;
    }
}

namespace lingdist
{
    DataFrame edit_dist_prepared(const PreparedInput &input, const IndexedCosts &cost, const FormPolicy &policy,
                                 bool detailed, const String &normalize_method, bool squareform,
                                 bool symmetric, int n_threads, bool quiet, bool cache)
    {
        validate_normalization(normalize_method);
        if (!input.symbol_table().same_id_space(cost.symbol_table()))
            throw std::invalid_argument("Input and costs must share the same SymbolTable.");
        const std::string method = normalize_method;
        const auto pairs = make_row_pairs(input, symmetric);
        StringVector labels1(pairs.size()), labels2(pairs.size());
        for (std::size_t i = 0; i < pairs.size(); ++i)
        {
            labels1[i] = Rcpp::String(input.row_names()[pairs[i].source], CE_UTF8);
            labels2[i] = Rcpp::String(input.row_names()[pairs[i].target], CE_UTF8);
        }

        // How to compute the distance for a single directed word pair from scratch (i.e. not from cache).
        const auto compute = [&](WordId a, WordId b)
        { return ed_dp(input.word(a), input.word(b), cost, method); };

        // We put all the distance computation logic into `distance_fn` and in this function, we concentrate on how to aggregate the results to get a final distance matrix.
        const auto aggregate = [&](const auto &distance_fn)
        {
            // Mean-only output needs one accumulator per row pair, not an O(row_pairs * columns) table of distances that will be discarded.
            std::vector<std::vector<double>> distances(detailed ? pairs.size() : 0);
            std::vector<double> averages(detailed ? 0 : pairs.size(), NA_REAL);
            std::unique_ptr<SafeProgressBar> bar;
            if (!quiet)
                bar = std::make_unique<SafeProgressBar>(pairs.size(), 1);
            parallel_for_checked(pairs.size(), [&](std::size_t i)
                                 {
                if (detailed) distances[i].resize(input.ncol());
                double sum = 0.0;
                std::size_t count = 0;
                for (std::size_t col = 0; col < input.ncol(); ++col)
                {
                    const double value = evaluate_forms(input.cell_forms(pairs[i].source, col),
                        input.cell_forms(pairs[i].target, col), policy, distance_fn);
                    if (detailed)
                        distances[i][col] = value;
                    else if (!Rcpp::NumericVector::is_na(value))
                    {
                        sum += value;
                        ++count;
                    }
                }
                if (!detailed && count) averages[i] = sum / count;
                if (bar) (*bar)++; }, n_threads);
            if (detailed)
                return gen_dist_df(std::move(distances), labels1, labels2, true,
                                   input.column_names(), squareform, symmetric);
            DataFrame result = DataFrame::create(Named("lab1") = labels1, Named("lab2") = labels2,
                                                 Named("dist") = Rcpp::wrap(averages));
            return squareform ? lingdist::long2squareform(result, symmetric) : result;
        };

        if (cache)
        {
            const auto plan = build_distance_plan(input, pairs, policy);
            std::vector<double> values(plan.tasks.size());
            {
                std::unique_ptr<SafeProgressBar> bar;
                if (!quiet)
                {
                    Rprintf("Computing %llu unique directed word-pair distances.\n",
                            static_cast<unsigned long long>(plan.tasks.size()));
                    bar = std::make_unique<SafeProgressBar>(plan.tasks.size(), 1);
                }
                // Storage is fixed before workers start; each task writes one distinct slot.
                parallel_for_checked(plan.tasks.size(), [&](std::size_t i)
                                     {
                values[i] = compute(plan.tasks[i].source, plan.tasks[i].target);
                if (bar) (*bar)++; }, n_threads);
            }
            // Both the task map and values remain read-only throughout the second parallel pass.
            return aggregate([&](WordId a, WordId b)
                             { return values[plan.task_ids.at(make_word_pair_key(a, b))]; });
        }
        else
        {
            // No caching: compute each directed word pair from scratch, even if it is repeated.
            return aggregate(compute);
        }
    }

    DataFrame edit_dist_df(const DataFrame &data, Nullable<DataFrame> cost_mat, const String &delim,
                           bool detailed, const String &normalize_method, const String &form_strategy,
                           const String &form_delim, std::vector<double> form_weights,
                           bool squareform, bool symmetric, int n_threads, bool check_missing_cost,
                           bool quiet, double default_sub_cost, double default_ins_del_cost, bool cache)
    {
        validate_normalization(normalize_method);
        const auto strategy = parse_form_strategy(static_cast<std::string>(form_strategy));
        const auto input = prepare_input(data, delim, form_delim, strategy != FormStrategy::Off);
        const auto policy = make_form_policy(strategy, input.max_forms(), std::move(form_weights));
        const auto cost = cost_mat.isNull()
                              ? IndexedCosts::uniform(input.symbol_table(), default_sub_cost, default_ins_del_cost)
                              : read_costs(input.symbol_table(), DataFrame(cost_mat));
        if (check_missing_cost && cost_mat.isNotNull())
            warn_missing(input.symbol_table(), DataFrame(cost_mat));
        return edit_dist_prepared(input, cost, policy, detailed, normalize_method, squareform, symmetric, n_threads, quiet, cache);
    }
}
