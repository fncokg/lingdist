#include <RcppThread.h>
#include <math.h>

#include "helpers.hpp"
#include "cost_table.hpp"
#include "dp.hpp"
#include "alignment.hpp"
#include "dist.hpp"

using namespace Rcpp;

namespace
{

    template <typename T, typename S>
    inline void insert_or_increment(std::unordered_map<T, S> &map, const T &key, const S &step)
    {
        auto result = map.insert({key, step});
        // result.second is true if inserted, false if already existed
        if (!result.second)
        {
            result.first->second += step;
        }
    }

    struct PMIRowResult
    {
        std::vector<std::pair<std::string, std::string>> char_pairs;
        std::vector<double> dist;
    };

    // Two chars are unordered, i.e. (c1,c2) == (c2,c1).
    // This is because in the alignment, replacement of c1 with c2 is equivalent to replacement of c2 with c1, since we randomly choose one as source and the other as target.
    struct AlignmentCount
    {
        // We DO allow diagonal storage here but currently,
        // but, we DONOT consider identical character pairs in PMI calculation
        lingdist::StrVec syms;
        size_t nsyms;
        std::vector<double> pair_counts;
        std::vector<double> sym_counts;

        AlignmentCount(const lingdist::StrVec &syms) : syms(syms)
        {
            nsyms = syms.size();
            size_t size_lower_tri = (nsyms * (nsyms + 1)) / 2;
            // NOTE: initialize with 1.0 to avoid any zero counts
            pair_counts.resize(size_lower_tri, 1.0);
            sym_counts.resize(nsyms, 1.0);
        }

        void clear()
        {
            std::fill(pair_counts.begin(), pair_counts.end(), 1.0);
            std::fill(sym_counts.begin(), sym_counts.end(), 1.0);
        }

        inline size_t rc2idx(size_t r, size_t c) const
        {
            if (r < c)
                std::swap(r, c);
            // return index in lower-triangular matrix
            // There're r full rows before the r-th row, and the r-th row has (r + 1) elements. Therefore, we have r(r+1)/2 elements before this row:

            // $\sum_{i=0}^{r-1} (i + 1) = \frac{r(r+1)}{2}$

            return (r * (r + 1)) / 2 + c;
        }

        // inline std::pair<double, double> idx2rc(size_t idx) const
        // {
        //     // First, solve the max r where r(r+1)/2 <= idx
        //     // That is, r^2 + r - 2*idx <= 0
        //     // => r = floor(( -1 + sqrt(1 + 8*idx) ) / 2)
        //     // Then, c = idx - r(r+1)/2

        //     // since idx is postive, we can directly use type cast to size_t instead of std::floor
        //     size_t r = static_cast<size_t>((std::sqrt(8.0 * static_cast<double>(idx) + 1.0) - 1.0) / 2.0);
        //     size_t c = idx - (r * (r + 1)) / 2;
        //     return {r, c};
        // }

        inline void add_pair(const std::string &c1, const std::string &c2, double step)
        {
            auto it1 = std::find(syms.begin(), syms.end(), c1);
            auto it2 = std::find(syms.begin(), syms.end(), c2);
            if (it1 != syms.end() && it2 != syms.end())
            {
                size_t idx1 = static_cast<size_t>(std::distance(syms.begin(), it1));
                size_t idx2 = static_cast<size_t>(std::distance(syms.begin(), it2));
                pair_counts[rc2idx(idx1, idx2)] += step;
            }
        }

        inline void add_sym(const std::string &c, double step)
        {
            auto it = std::find(syms.begin(), syms.end(), c);
            if (it != syms.end())
            {
                size_t idx = static_cast<size_t>(std::distance(syms.begin(), it));
                sym_counts[idx] += step;
            }
        }
    };

    void update_cost_table(const std::vector<lingdist::AlignmentResult> &results, lingdist::CostTable &cost, AlignmentCount &alignment_count)
    {
        alignment_count.clear();
        // we do not clear cost, since all its entries will be updated anyway
        // double total_count = 0.0;
        double eff_val;

        for (const auto &res : results)
        {
            for (const auto &path : res.alignments)
            {
                if (path.empty())
                    continue;
                // This is not mentioned in Wieling(2012), but seems necessary:
                // If we find multiple optimal paths, we must make sure that pairs in each path only contribute fractionally,
                // so that this alignment contributes only 1.0 to the all counts.
                eff_val = 1.0 / static_cast<double>(path.size());
                for (const auto &char_pair : path)
                {

                    const std::string &c1 = char_pair.first;
                    const std::string &c2 = char_pair.second;
                    // We now only consider non-identical character pairs
                    if (c1 != c2)
                    {
                        alignment_count.add_pair(c1, c2, eff_val);
                        alignment_count.add_sym(c1, eff_val);
                        alignment_count.add_sym(c2, eff_val);
                        // total_count += eff_val;
                    }
                }
            }
        }

        // Now we loop through all the lower-triangular entries twice:
        // 1st time to get raw distances and find min/max
        // 2nd time to normalize and update cost table
        // note that in both loops, we skip diagonal entries (identical character pairs)

        double max_raw_dist = -std::numeric_limits<double>::infinity();
        double min_raw_dist = std::numeric_limits<double>::infinity();

        for (size_t r = 0; r < alignment_count.nsyms; r++)
        {
            // skip identical pairs
            for (size_t c = 0; c < r; c++)
            {
                size_t idx = alignment_count.rc2idx(r, c);
                double c_xy = alignment_count.pair_counts[idx];
                double c_x = alignment_count.sym_counts[r];
                double c_y = alignment_count.sym_counts[c];

                // $N$: total number of all character pairs aligned
                // $c(x,y)$: count of co-occurrence of x and y
                // $c(x),c(y)$: count of occurrence of x and y respectively

                // $\mathrm{PMI} = \log_2 \frac{p(x,y)}{p(x)p(y)}=\log_2 \frac{c(x,y)N}{c(x)c(y)}$

                // Thus $\mathrm{PMI} = \frac1{\ln2}[\ln \frac{c(x,y)}{c(x)c(y)}+ \ln N]$

                // $N$ and $\ln 2$ are constants participating only in a linear transformation, which means
                // after normalization, they do not affect the relative distances.
                // we only need to update the cost table with - log(c(x, y) / (c(x) c(y)))

                // Directly log is safe here since we have initialized counts with 1.0 to avoid zero counts.
                double raw_dist = -std::log((c_xy) / (c_x * c_y));
                max_raw_dist = std::max(max_raw_dist, raw_dist);
                min_raw_dist = std::min(min_raw_dist, raw_dist);

                // reuse pair_counts to store raw_dist
                alignment_count.pair_counts[idx] = raw_dist;
            }
        }

        for (size_t r = 0; r < alignment_count.nsyms; r++)
        {
            // we do not skip diagonal entries
            // we explicitly set them to 0.0 cost later
            for (size_t c = 0; c <= r; c++)
            {
                const std::string &c1 = alignment_count.syms[r];
                const std::string &c2 = alignment_count.syms[c];
                if (r == c)
                {
                    cost.set_cost(c1, c2, 0.0);
                    continue;
                }
                size_t idx = alignment_count.rc2idx(r, c);
                double raw_dist = alignment_count.pair_counts[idx];
                double norm_dist;
                if (lingdist::double_equal(max_raw_dist, min_raw_dist))
                {

                    norm_dist = 1.0; // all distances are identical
                }
                else
                {
                    norm_dist = (raw_dist - min_raw_dist) / (max_raw_dist - min_raw_dist);
                }
                cost.set_cost(c1, c2, norm_dist);
                cost.set_cost(c2, c1, norm_dist);
            }
        }
    }

    lingdist::AlignmentResult pmi_row(const std::vector<lingdist::StrVec> &row1, const std::vector<lingdist::StrVec> &row2, const lingdist::CostTable &cost, int alignment_max_paths = 3)
    {
        lingdist::AlignmentResult result;
        std::size_t ncols = row1.size();
        double sum_dist = 0.0, nwords = 0.0;
        lingdist::AlignmentResult this_res;
        for (std::size_t coli = 0; coli < ncols; coli++)
        {
            const auto &chars1 = row1[coli];
            const auto &chars2 = row2[coli];
            if (!chars1.empty() && !chars2.empty())
            {
                this_res = lingdist::get_string_alignment_result(chars1, chars2, cost, alignment_max_paths);
                sum_dist += this_res.distance;
                nwords += 1.0;
                result.alignments.insert(result.alignments.end(),
                                         this_res.alignments.begin(),
                                         this_res.alignments.end());
            }
        }
        result.distance = nwords > 0 ? sum_dist / nwords : NA_REAL;
        return result;
    }
}

namespace lingdist
{
    List pmi_df(const DataFrame &data, const String &delim, bool squareform, bool parallel, int n_threads, int max_epochs, double tol, int alignment_max_paths, bool quiet)
    {
        if (!quiet)
            Rprintf("Starting PMI distance computation on data frame with %d rows...\n", data.nrow());

        // First, prepare a default cost table
        lingdist::StrVec unique_chars = lingdist::get_all_unique_syms(data, delim, true);

        if (!quiet)
            Rprintf("Identified %d unique characters in the data frame.\n", static_cast<int>(unique_chars.size()));
        lingdist::CostTable prev_cost, cost = lingdist::build_default_cost_table(unique_chars);
        if (!quiet)
            Rprintf("Initialized default cost table.\n");

        List report;
        AlignmentCount alignment_count(unique_chars);
        auto rows_vector = lingdist::split_df(data, delim);
        auto [lab1Col, lab2Col, row_pairs] = lingdist::get_row_pairs(data, true);

        int n_row_pairs = static_cast<int>(row_pairs.size());
        if (!quiet)
        {
            Rprintf("Starting PMI iterations with %d row pairs...\n", n_row_pairs);
            if (parallel)
                Rprintf("Using %d threads for parallel computation.\n", n_threads);
            else
                Rprintf("Using single-threaded computation.\n");
        }
        double sum_diff = 0.0, mean_diff = 0.0;
        std::vector<lingdist::AlignmentResult> results(n_row_pairs);
        int iepoch;
        for (iepoch = 1; iepoch <= max_epochs; iepoch++)
        {
            if (!quiet)
                Rprintf("Epoch %d/%d\n", iepoch, max_epochs);

            // RcppThread::ProgressBar *bar = nullptr;
            std::unique_ptr<lingdist::SafeProgressBar> bar;
            if (!quiet)
                bar = std::make_unique<lingdist::SafeProgressBar>(n_row_pairs, 1);
            std::function<void(int)> loop_body = [&](int idx)
            {
                auto [rowi, rowj] = row_pairs[idx];
                results[idx] = pmi_row(rows_vector[rowi], rows_vector[rowj], cost, alignment_max_paths);
                if (bar)
                    (*bar)++;
            };
            if (parallel)
            {
                RcppThread::parallelFor(0, n_row_pairs, loop_body, n_threads);
            }
            else
            {
                lingdist::singleFor(0, n_row_pairs, loop_body);
            }

            if (!quiet)
                Rprintf("Updating cost table...\n");
            prev_cost = cost;
            update_cost_table(results, cost, alignment_count);
            // Check convergence
            sum_diff = 0.0;
            for (size_t idx = 0; idx < cost.data.size(); ++idx)
            {
                sum_diff += std::fabs(cost.data[idx] - prev_cost.data[idx]);
            }
            mean_diff = sum_diff / (static_cast<double>(cost.data.size()) - static_cast<double>(unique_chars.size())); // exclude diagonal entries
            if (!quiet)
                Rprintf("Cost table updated, mean absolute difference: %.6f, total absolute difference: %.6f\n", mean_diff, sum_diff);
            if (mean_diff < tol)
            {
                if (!quiet)
                    Rprintf("Convergence reached with mean difference %.6f < tol %.6f. Stopping iterations.\n", mean_diff, tol);
                break;
            }
        }
        if (mean_diff >= tol)
        {
            warning("PMI distance computation did not converge within the maximum number of epochs. Try increasing max_epochs or tol.");
        }

        std::vector<double> dists(row_pairs.size());
        for (std::size_t idx = 0; idx < row_pairs.size(); idx++)
        {
            dists[idx] = results[idx].distance;
        }
        DataFrame result = DataFrame::create(Named("lab1") = lab1Col, Named("lab2") = lab2Col, Named("dist") = NumericVector::import(dists.begin(), dists.end()));
        if (squareform)
        {
            result = lingdist::long2squareform(result, true);
        }
        report["result"] = result;
        // Note: when we finish the iteration, i.e. we finish the last update of cost table, the cost table we used to compute distances is actually the previous one.
        report["cost"] = prev_cost.to_dataframe();
        report["sum_diff"] = sum_diff;
        report["mean_diff"] = mean_diff;
        report["converged"] = mean_diff < tol;
        report["n_epochs"] = iepoch;
        if (!quiet)
            Rprintf("PMI distance computation completed.\n");
        return report;
    }
}