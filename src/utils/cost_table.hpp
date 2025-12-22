#pragma once

#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include "helpers.hpp"

using namespace Rcpp;

namespace lingdist
{

    // 高性能成本表：列名/行名分别映射到索引，数据按列主序或行主序存放（此处采用 col-major: idx = r * ncol + c）
    struct CostTable
    {
        std::unordered_map<std::string, int> colIndex; // 符号 -> 列索引（对应被删除/插入一侧）
        std::unordered_map<std::string, int> rowIndex; // 符号 -> 行索引
        std::vector<double> data;                      // 长度 nrow * ncol
        int nrow{0};
        int ncol{0};

        double get_cost(const std::string &str1, const std::string &str2) const;

        void set_cost(const std::string &str1, const std::string &str2, double value);

        DataFrame to_dataframe() const;

        inline bool has(const std::string &cname, const std::string &rname) const
        {
            return colIndex.find(cname) != colIndex.end() && rowIndex.find(rname) != rowIndex.end();
        }
        inline size_t get_index(int c, int r) const
        {
            return static_cast<size_t>(r) * static_cast<size_t>(ncol) + static_cast<size_t>(c);
        }

        inline double at_by_index(int c, int r) const
        {
            return data[get_index(c, r)];
        }
    };

    // 从 DataFrame 构建 CostTable（允许列名与行名集合不同）
    CostTable build_cost_table(const Rcpp::DataFrame &cost_mat);

    CostTable build_default_cost_table(const lingdist::StrVec &chars, double sub_cost = 1.0, double ins_del_cost = 1.0);

} // namespace lingdist
