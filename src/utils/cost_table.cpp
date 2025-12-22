#include <Rcpp.h>

#include "cost_table.hpp"

using namespace Rcpp;

namespace lingdist
{

    CostTable build_cost_table(const DataFrame &cost_mat)
    {
        CostTable table;
        if (cost_mat.nrow() == 0 || cost_mat.ncol() == 0)
        {
            return table;
        }
        StringVector _row_names = cost_mat.attr("row.names");
        StringVector _col_names = cost_mat.names();
        StrVec row_names = as<StrVec>(_row_names);
        StrVec col_names = as<StrVec>(_col_names);

        table.nrow = static_cast<int>(row_names.size());
        table.ncol = static_cast<int>(col_names.size());
        table.data.resize(static_cast<size_t>(table.nrow) * static_cast<size_t>(table.ncol));

        for (int c = 0; c < table.ncol; ++c)
            table.colIndex[col_names[c]] = c;
        for (int r = 0; r < table.nrow; ++r)
            table.rowIndex[row_names[r]] = r;

        for (int c = 0; c < table.ncol; ++c)
        {
            NumericVector this_col = cost_mat[c];
            for (int r = 0; r < table.nrow; ++r)
            {
                table.data[table.get_index(c, r)] = this_col[r];
            }
        }
        return table;
    }

    CostTable build_default_cost_table(const lingdist::StrVec &chars, double sub_cost, double ins_del_cost)
    {
        CostTable table;
        int nchars = static_cast<int>(chars.size());
        table.nrow = nchars;
        table.ncol = nchars;
        table.data.resize(static_cast<size_t>(table.nrow) * static_cast<size_t>(table.ncol));
        for (int i = 0; i < nchars; ++i)
        {
            table.rowIndex[chars[i]] = i;
            table.colIndex[chars[i]] = i;
        }
        for (int r = 0; r < table.nrow; ++r)
        {
            for (int c = 0; c < table.ncol; ++c)
            {
                if (r == c)
                {
                    table.data[table.get_index(c, r)] = 0.0;
                }
                else if (r == 0 || c == 0)
                {
                    table.data[table.get_index(c, r)] = ins_del_cost; // insertion/deletion cost
                }
                else
                {
                    table.data[table.get_index(c, r)] = sub_cost;
                }
            }
        }
        return table;
    }

    double CostTable::get_cost(const std::string &str1, const std::string &str2) const
    {
        auto cIt = colIndex.find(str1);
        auto rIt = rowIndex.find(str2);
        if (cIt != colIndex.end() && rIt != rowIndex.end())
        {
            return at_by_index(cIt->second, rIt->second);
        }
        // 未定义时，遵循原逻辑：相同为 0，否则 1
        return (str1 == str2) ? 0.0 : 1.0;
    }

    void CostTable::set_cost(const std::string &str1, const std::string &str2, double new_value)
    {
        auto cIt = colIndex.find(str1);
        auto rIt = rowIndex.find(str2);
        if (cIt != colIndex.end() && rIt != rowIndex.end())
        {
            data[get_index(cIt->second, rIt->second)] = new_value;
        }
    }

    DataFrame CostTable::to_dataframe() const
    {
        StringVector row_names(nrow);
        StringVector col_names(ncol);
        for (const auto &pair : rowIndex)
        {
            row_names[pair.second] = pair.first;
        }
        for (const auto &pair : colIndex)
        {
            col_names[pair.second] = pair.first;
        }
        DataFrame df = DataFrame::create();
        for (int c = 0; c < ncol; ++c)
        {
            NumericVector this_col(nrow);
            for (int r = 0; r < nrow; ++r)
            {
                this_col[r] = at_by_index(c, r);
            }
            df.push_back(this_col);
        }
        df.attr("names") = col_names;
        df.attr("row.names") = row_names;
        return df;
    }

} // namespace lingdist
