#include <Rcpp.h>

#include "cost_table.hpp"

using namespace Rcpp;

namespace lingdist
{
    lingdist::StrVec CostTable::check_missing_symbols(const DataFrame &data, const String &delim) const
    {
        lingdist::StrVec unique_chars = lingdist::get_all_unique_chars(data, delim);
        unique_chars.push_back(lingdist::EMPTY);
        lingdist::StrVec missing_chars;
        for (const auto &ch : unique_chars)
        {
            if (sym_index.find(ch) == sym_index.end())
            {
                missing_chars.push_back(ch);
            }
        }
        return missing_chars;
    }

    CostTable build_cost_table(const DataFrame &cost_mat)
    {
        CostTable table;
        if (cost_mat.nrow() == 0 || cost_mat.ncol() == 0)
        {
            return table;
        }
        if (cost_mat.nrow() != cost_mat.ncol())
        {
            stop("Cost matrix must be square.");
        }

        StrVec row_names = as<StrVec>(cost_mat.attr("row.names"));
        StrVec col_names = as<StrVec>(cost_mat.names());

        table.is_fast = false;
        table.nsyms = static_cast<int>(row_names.size());
        table.syms = row_names;
        table.data.resize(static_cast<size_t>(table.nsyms) * static_cast<size_t>(table.nsyms));

        std::unordered_map<std::string, int> cindexer, rindexer;
        for (int i = 0; i < table.nsyms; ++i)
        {
            cindexer[col_names[i]] = i;
            rindexer[row_names[i]] = i;
        }
        table.sym_index = rindexer;

        for (auto &sym1 : table.syms)
        {
            int c_id = cindexer[sym1];
            NumericVector this_col = cost_mat[c_id];
            for (auto &sym2 : table.syms)
            {
                int r_id = rindexer[sym2];
                table.set_cost(sym1, sym2, this_col[r_id]);
            }
        }
        return table;
    }

    CostTable build_fast_cost_table(const lingdist::StrVec &syms, double sub_cost, double ins_del_cost)
    {
        CostTable table;
        table.is_fast = true;
        table.fast_sub_cost = sub_cost;
        table.fast_ins_del_cost = ins_del_cost;
        // copy
        table.syms = syms;
        return table;
    }

    CostTable build_default_cost_table(const lingdist::StrVec &syms, double sub_cost, double ins_del_cost)
    {
        CostTable table;
        int nsyms = static_cast<int>(syms.size());
        table.is_fast = false;
        table.syms = syms;
        table.nsyms = nsyms;
        table.data.resize(static_cast<size_t>(table.nsyms) * static_cast<size_t>(table.nsyms));
        for (int i = 0; i < nsyms; ++i)
        {
            table.sym_index[syms[i]] = i;
        }
        for (auto &sym1 : syms)
        {
            for (auto &sym2 : syms)
            {
                double value = (sym1 == sym2) ? 0.0 : sub_cost;
                if (sym1 == lingdist::EMPTY || sym2 == lingdist::EMPTY)
                {
                    value = ins_del_cost;
                }
                table.set_cost(sym1, sym2, value);
            }
        }
        return table;
    }

    double CostTable::get_cost(const std::string &str1, const std::string &str2) const
    {
        if (is_fast)
        {
            if (str1 == str2)
            {
                return 0.0;
            }
            else if (str1 == lingdist::EMPTY || str2 == lingdist::EMPTY)
            {
                return fast_ins_del_cost;
            }
            else
            {
                return fast_sub_cost;
            }
        }
        else
        {
            auto cfind = sym_index.find(str1);
            auto rfind = sym_index.find(str2);
            if (cfind != sym_index.end() && rfind != sym_index.end())
            {
                return at_by_index(cfind->second, rfind->second);
            }
            // 未定义时，遵循原逻辑：相同为 0，否则 1
            return (str1 == str2) ? 0.0 : 1.0;
        }
    }

    void CostTable::set_cost(const std::string &str1, const std::string &str2, double new_value)
    {
        auto cfind = sym_index.find(str1);
        auto rfind = sym_index.find(str2);
        if (cfind != sym_index.end() && rfind != sym_index.end())
        {
            data[get_index(cfind->second, rfind->second)] = new_value;
        }
    }

    DataFrame CostTable::to_dataframe() const
    {
        if (is_fast)
        {
            return build_default_cost_table(
                       syms, fast_sub_cost, fast_ins_del_cost)
                .to_dataframe();
        }
        else
        {
            StringVector row_names(nsyms);
            StringVector col_names(nsyms);
            for (const auto &[key, value] : sym_index)
            {
                row_names[value] = key;
                col_names[value] = key;
            }
            DataFrame df = DataFrame::create();
            for (int c = 0; c < nsyms; ++c)
            {
                NumericVector this_col(nsyms);
                for (int r = 0; r < nsyms; ++r)
                {
                    this_col[r] = at_by_index(c, r);
                }
                df.push_back(this_col);
            }
            df.attr("names") = col_names;
            df.attr("row.names") = row_names;
            return df;
        }
    }

} // namespace lingdist
