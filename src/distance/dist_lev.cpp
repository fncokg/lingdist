#include <Rcpp.h>
#include <RcppThread.h>

#include "dist.hpp"
#include "dp.hpp"
#include "dispatcher.hpp"
#include "cost_io.hpp"
#include "prepare_input.hpp"

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
}

namespace lingdist
{
    // TODO: Utilize aggregation or caching to avoid repeated evaluation of the same form pairs.
    DataFrame edit_dist_prepared(const PreparedInput &input, const IndexedCosts &cost, const FormPolicy &policy,
                                 bool detailed, const String &normalize_method, bool squareform,
                                 bool symmetric, int n_threads, bool quiet)
    {
        validate_normalization(normalize_method);
        if (!input.symbol_table().same_id_space(cost.symbol_table()))
            throw std::invalid_argument("Input and costs must share the same SymbolTable.");
        const std::string method = normalize_method;
        const auto pairs = make_row_pairs(input, symmetric);
        StringVector labels1(pairs.size()), labels2(pairs.size());
        for (std::size_t i = 0; i < pairs.size(); ++i)
        {
            labels1[i] = input.row_names()[pairs[i].source];
            labels2[i] = input.row_names()[pairs[i].target];
        }
        std::vector<std::vector<double>> distances(pairs.size());
        std::unique_ptr<SafeProgressBar> bar;
        if (!quiet)
            bar = std::make_unique<SafeProgressBar>(pairs.size(), 1);
        RcppThread::parallelFor(0, pairs.size(), [&](std::size_t i)
                                {
            auto &row = distances[i];
            row.assign(input.ncol(), NA_REAL);
            for (std::size_t col = 0; col < input.ncol(); ++col)
            {
                row[col] = evaluate_forms(input.cell_forms(pairs[i].source, col),
                    input.cell_forms(pairs[i].target, col), policy, [&](WordId a, WordId b)
                    { return ed_dp(input.word(a), input.word(b), cost, method); });
            }
            if (bar) (*bar)++; }, n_threads);
        return gen_dist_df(std::move(distances), labels1, labels2, detailed,
                           input.column_names(), squareform, symmetric);
    }

    DataFrame edit_dist_df(const DataFrame &data, Nullable<DataFrame> cost_mat, const String &delim,
                           bool detailed, const String &normalize_method, const String &form_strategy,
                           const String &form_delim, std::vector<double> form_weights,
                           bool squareform, bool symmetric, int n_threads, bool check_missing_cost,
                           bool quiet, double default_sub_cost, double default_ins_del_cost)
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
        return edit_dist_prepared(input, cost, policy, detailed, normalize_method, squareform, symmetric, n_threads, quiet);
    }
}
