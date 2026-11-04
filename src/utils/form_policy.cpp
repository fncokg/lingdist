#include "form_policy.hpp"

#include <cmath>
#include <utility>

namespace lingdist
{
    FormStrategy parse_form_strategy(const std::string &name)
    {
        if (name == "off") return FormStrategy::Off;
        if (name == "first") return FormStrategy::First;
        if (name == "weighting") return FormStrategy::Weighting;
        if (name == "mean_all") return FormStrategy::MeanAll;
        if (name == "min_all") return FormStrategy::Minimum;
        throw std::invalid_argument("Unsupported form_strategy. Use off, first, weighting, mean_all or min_all.");
    }

    FormPolicy make_form_policy(FormStrategy strategy, std::size_t max_forms, std::vector<double> weights)
    {
        if (strategy != FormStrategy::Weighting) return {strategy, {}};
        if (weights.empty())
        {
            weights.resize(max_forms);
            double weight = 1.0;
            for (auto &value : weights) { value = weight; weight *= 0.5; }
        }
        if (weights.size() < max_forms)
            throw std::invalid_argument("form_weights must cover all observed forms.");
        for (const auto weight : weights)
            if (!std::isfinite(weight) || weight < 0.0)
                throw std::invalid_argument("form_weights must be finite and non-negative.");
        return {strategy, std::move(weights)};
    }
}
