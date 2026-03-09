#include "helpers.hpp"

double dispatcher_weighting(const std::vector<lingdist::StrVec> &forms1, const std::vector<lingdist::StrVec> &forms2, std::function<double(const lingdist::StrVec &, const lingdist::StrVec &)> dist_func, const std::vector<double> &weights);

double dispatcher_all(const std::vector<lingdist::StrVec> &forms1, const std::vector<lingdist::StrVec> &forms2, std::function<double(const lingdist::StrVec &, const lingdist::StrVec &)> dist_func, int mode);