#pragma once

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace lingdist
{
    enum class FormStrategy { Off, First, Weighting, MeanAll, Minimum };

    struct FormPolicy
    {
        FormStrategy strategy;
        std::vector<double> weights;
    };

    FormStrategy parse_form_strategy(const std::string &name);
    FormPolicy make_form_policy(FormStrategy strategy, std::size_t max_forms,
                                std::vector<double> weights = {});

    // Shared by training planning and prediction. Every valid cell pair emits
    // total mass 1, BEFORE aggregation across observations. Empty cells and
    // zero effective weight are invalid. Padding repeats a form, never a gap.
    // Containers may be vectors of WJD sequences or PreparedInput form views.
    template <class FormsA, class FormsB, class Emit>
    bool for_each_form_pair(const FormsA &a, const FormsB &b,
                            const FormPolicy &policy, Emit emit)
    {
        if (a.empty() || b.empty()) return false;
        if (policy.strategy == FormStrategy::Off || policy.strategy == FormStrategy::First)
        {
            emit(a[0], b[0], 1.0);
            return true;
        }
        if (policy.strategy == FormStrategy::Weighting)
        {
            const auto length = std::max(a.size(), b.size());
            if (policy.weights.size() < length)
                throw std::invalid_argument("form_weights must cover all observed forms.");
            // Scaling avoids overflow even for valid weights near DBL_MAX.
            const double scale = *std::max_element(policy.weights.begin(), policy.weights.begin() + length);
            if (scale == 0.0) return false;
            double total = 0.0;
            for (std::size_t i = 0; i < length; ++i) total += policy.weights[i] / scale;
            for (std::size_t i = 0; i < length; ++i)
                if (policy.weights[i] > 0.0)
                {
                    const double weight = (policy.weights[i] / scale) / total;
                    // A normalized subnormal can underflow; never materialize
                    // zero-weight tasks or count them as emitted form pairs.
                    if (weight > 0.0)
                        emit(a[std::min(i, a.size() - 1)], b[std::min(i, b.size() - 1)], weight);
                }
            return true;
        }
        // Minimum has the same candidates as MeanAll, but prediction takes their
        // minimum. It cannot form a fixed, cost-independent PMI training plan.
        const double weight = 1.0 / static_cast<double>(a.size()) / static_cast<double>(b.size());
        for (std::size_t i = 0; i < a.size(); ++i)
            for (std::size_t j = 0; j < b.size(); ++j)
                emit(a[i], b[j], weight);
        return true;
    }
}
