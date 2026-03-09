#include "helpers.hpp"

double dispatcher_weighting(const std::vector<lingdist::StrVec> &forms1, const std::vector<lingdist::StrVec> &forms2, lingdist::DistFunc dist_func, const std::vector<double> &weights);

double dispatcher_all(const std::vector<lingdist::StrVec> &forms1, const std::vector<lingdist::StrVec> &forms2, lingdist::DistFunc dist_func, int mode);

double dispatcher_first(const std::vector<lingdist::StrVec> &forms1, const std::vector<lingdist::StrVec> &forms2, lingdist::DistFunc dist_func);

lingdist::FormsDispatcher make_dispatcher(const String &forms_strategy, lingdist::DistFunc dist_func, const std::vector<double> &forms_weights);