#pragma once

#include <cmath>

namespace lingdist
{
    constexpr double negative_infinity = -std::numeric_limits<double>::infinity();

    inline constexpr double EPS = 1e-9;

    // This is a local DP tie rule, not a tolerance on complete path costs.
    inline bool double_equal(double a, double b) noexcept
    {
        return std::abs(a - b) < EPS;
    }

    /*
    $\mathrm{logaddexp}(a,b)=c \Rightarrow \log(\exp(a)+\exp(b))=c$
    We use this to get the sum of two log probabilities without underflowing to zero.
    The result is always finite if at least one of the inputs is finite.
    */
    inline double logaddexp(double a, double b)
    {
        if (a == negative_infinity)
            return b;
        if (b == negative_infinity)
            return a;
        const auto larger = std::max(a, b);
        return larger + std::log1p(std::exp(-std::abs(a - b)));
    }
}
