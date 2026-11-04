#pragma once

#include "helpers.hpp"
#include "form_policy.hpp"

// The same normalized pair generator feeds the training histogram and distance
// reducers. Only min_all uses a non-linear reducer and is prediction-only.
template <class FormsA, class FormsB, class Distance>
double evaluate_forms(const FormsA &a, const FormsB &b,
                      const lingdist::FormPolicy &policy, Distance distance)
{
    double sum = 0.0, mass = 0.0;
    double minimum = std::numeric_limits<double>::infinity();
    const bool skip_missing = policy.strategy == lingdist::FormStrategy::MeanAll ||
                              policy.strategy == lingdist::FormStrategy::Minimum;
    lingdist::for_each_form_pair(a, b, policy, [&](const auto &x, const auto &y, double weight)
    {
        const double value = distance(x, y);
        if (skip_missing && Rcpp::NumericVector::is_na(value)) return;
        sum += weight * value;
        mass += weight;
        minimum = std::min(minimum, value);
    });
    if (mass == 0.0) return NA_REAL;
    return policy.strategy == lingdist::FormStrategy::Minimum ? minimum : sum / mass;
}

template <class Sequence>
using SequenceDistance = std::function<double(const Sequence &, const Sequence &)>;
template <class Sequence>
using SequenceDispatcher = std::function<double(const std::vector<Sequence> &, const std::vector<Sequence> &)>;

template <class Sequence>
SequenceDispatcher<Sequence> make_dispatcher(const String &strategy, SequenceDistance<Sequence> distance,
                                            const std::vector<double> &weights)
{
    const auto policy = lingdist::make_form_policy(
        lingdist::parse_form_strategy(static_cast<std::string>(strategy)), weights.size(), weights);
    return [policy, distance](const auto &a, const auto &b)
    { return evaluate_forms(a, b, policy, distance); };
}
