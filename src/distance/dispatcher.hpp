#include <vector>
#include <functional>

template <typename T>
double dispatcher_weighting(const std::vector<T> &forms1, const std::vector<T> &forms2, std::function<double(const T &, const T &)> dist_func, std::vector<double> weights);

template <typename T>
double dispatcher_all(const std::vector<T> &forms1, const std::vector<T> &forms2, std::function<double(const T &, const T &)> dist_func, int mode);