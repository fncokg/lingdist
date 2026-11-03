#pragma once

#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include "helpers.hpp"

using namespace Rcpp;

namespace lingdist
{

    // All accessors take (source, target); storage index = target * nsyms + source.
    struct CostTable
    {
        // A table structure
        std::unordered_map<std::string, int> sym_index;
        lingdist::StrVec syms;
        int nsyms{0};
        std::vector<double> data; // 长度 nsyms * nsyms

        // A fast table structure
        bool is_fast{false}; // 是否为快速默认成本表
        double fast_sub_cost{1.0};
        double fast_ins_del_cost{1.0};

        double get_cost(const std::string &str1, const std::string &str2) const;

        void set_cost(const std::string &str1, const std::string &str2, double value);

        DataFrame to_dataframe() const;

        inline bool has(const std::string &cname, const std::string &rname) const
        {
            return sym_index.find(cname) != sym_index.end() && sym_index.find(rname) != sym_index.end();
        }
        inline size_t get_index(int source, int target) const
        {
            return static_cast<size_t>(target) * static_cast<size_t>(nsyms) + static_cast<size_t>(source);
        }

        inline double at_by_index(int source, int target) const
        {
            return data[get_index(source, target)];
        }

        lingdist::StrVec check_missing_symbols(const DataFrame &data, const String &delim) const;
    };

    // Build a cost table from a DataFrame with matching row/column symbol sets.
    CostTable build_cost_table(const Rcpp::DataFrame &cost_mat);

    CostTable build_fast_cost_table(const lingdist::StrVec &syms, double sub_cost = 1.0, double ins_del_cost = 1.0);

    CostTable build_default_cost_table(const lingdist::StrVec &syms, double sub_cost = 1.0, double ins_del_cost = 1.0);

} // namespace lingdist
