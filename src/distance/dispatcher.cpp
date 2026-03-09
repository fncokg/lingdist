#include "dispatcher.hpp"

double dispatcher_weighting(const std::vector<lingdist::StrVec> &forms1, const std::vector<lingdist::StrVec> &forms2, lingdist::DistFunc dist_func, const std::vector<double> &weights)
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
double dispatcher_all(const std::vector<lingdist::StrVec> &forms1, const std::vector<lingdist::StrVec> &forms2, lingdist::DistFunc dist_func, int mode)
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

double dispatcher_first(const std::vector<lingdist::StrVec> &forms1, const std::vector<lingdist::StrVec> &forms2, lingdist::DistFunc dist_func)
{
    if (forms1.empty() || forms2.empty())
        return NA_REAL;
    return dist_func(forms1[0], forms2[0]);
}

lingdist::FormsDispatcher make_dispatcher(const String &forms_strategy, lingdist::DistFunc dist_func, const std::vector<double> &forms_weights)
{
    if (forms_strategy == "weighting")
    {
        return [dist_func, forms_weights](const std::vector<lingdist::StrVec> &forms1, const std::vector<lingdist::StrVec> &forms2)
        { return dispatcher_weighting(forms1, forms2, dist_func, forms_weights); };
    }
    else if (forms_strategy == "mean_all")
    {
        return [dist_func](const std::vector<lingdist::StrVec> &forms1, const std::vector<lingdist::StrVec> &forms2)
        { return dispatcher_all(forms1, forms2, dist_func, 0); };
    }
    else if (forms_strategy == "min_all")
    {
        return [dist_func](const std::vector<lingdist::StrVec> &forms1, const std::vector<lingdist::StrVec> &forms2)
        { return dispatcher_all(forms1, forms2, dist_func, 1); };
    }
    else if (forms_strategy == "off" || forms_strategy == "first")
    {
        return [dist_func](const std::vector<lingdist::StrVec> &forms1, const std::vector<lingdist::StrVec> &forms2)
        { return dispatcher_first(forms1, forms2, dist_func); };
    }
    else
    {
        stop("Unsupported forms_strategy. Supported values are 'weighting', 'mean_all', 'min_all', 'off' and 'first'.");
    }
}