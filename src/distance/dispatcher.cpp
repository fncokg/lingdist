#include "dispatcher.hpp"

template <typename T>
double dispatcher_weighting(const std::vector<T> &forms1, const std::vector<T> &forms2, std::function<double(const T &, const T &)> dist_func, std::vector<double> weights)
{
    double sum_dist = 0.0, sum_weights = 0.0;
    size_t nforms1 = forms1.size(), nforms2 = forms2.size();
    if (nforms1 != 0 && nforms2 != 0)
    {
        size_t max_len = std::max(nforms1, nforms2);
        for (size_t i = 0; i < max_len; i++)
        {
            const auto &form1 = i < nforms1 ? forms1[i] : forms1.back();
            const auto &form2 = i < nforms2 ? forms2[i] : forms2.back();
            sum_dist += dist_func(form1, form2) * weights[i];
            sum_weights += weights[i];
        }
    }
    if (sum_weights <= 0.0)
        return NA_REAL;
    return sum_dist / sum_weights;
}

// mode: 0 for mean, 1 for min
template <typename T>
double dispatcher_all(const std::vector<T> &forms1, const std::vector<T> &forms2, std::function<double(const T &, const T &)> dist_func, int mode)
{
    std::vector<double> dists;
    size_t nforms1 = forms1.size(), nforms2 = forms2.size();
    if (nforms1 != 0 && nforms2 != 0)
    {
        for (size_t i = 0; i < nforms1; i++)
        {
            for (size_t j = 0; j < nforms2; j++)
            {
                double this_dist = dist_func(forms1[i], forms2[j]);
                if (!Rcpp::NumericVector::is_na(this_dist))
                {
                    dists.push_back(this_dist);
                }
            }
        }
    }
    if (dists.empty())
        return NA_REAL;
    if (mode == 1)
    {
        return *std::min_element(dists.begin(), dists.end());
    }
    else
    {
        return std::accumulate(dists.begin(), dists.end(), 0.0) / dists.size();
    }
}