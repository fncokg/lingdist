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

    // We need a structure to hold character pair counts
    // The key is that: two chars are unordered, i.e. (c1,c2) == (c2,c1).
    // This is because in the alignment, replacement of c1 with c2 is equivalent to replacement of c2 with c1, since we randomly choose one as source and the other as target.
    struct CharMap
    {
        std::unordered_map<std::string, uint32_t> idx;
        std::unordered_map<uint32_t, std::string> rev_idx;

        std::unordered_map<uint64_t, double> pair_count;

        CharMap(const std::vector<std::string> &chars)
        {
            uint32_t index = 0;
            for (const auto &c : chars)
            {
                idx[c] = index;
                rev_idx[index] = c;
                index++;
            }
        }

        inline bool has(std::string c1, std::string c2) const
        {
            return pair_count.find(s2k(c1, c2)) != pair_count.end();
        }

        inline std::pair<std::string, std::string> k2s(uint64_t key) const
        {
            uint32_t i1 = static_cast<uint32_t>(key >> 32);
            uint32_t i2 = static_cast<uint32_t>(key & 0xFFFFFFFF);
            return {rev_idx.at(i1), rev_idx.at(i2)};
        }

        inline uint64_t s2k(std::string c1, std::string c2) const
        {
            // Unsafe if c1 or c2 not in idx
            uint32_t i1 = idx.at(c1);
            uint32_t i2 = idx.at(c2);
            // let smaller index be in higher bits to ensure uniqueness
            return i1 <= i2 ? (static_cast<uint64_t>(i1) << 32) | i2 : (static_cast<uint64_t>(i2) << 32) | i1;
        }

        inline void increment_item(std::string c1, std::string c2, double step)
        {
            insert_or_increment(pair_count, s2k(c1, c2), step);
        }
    };

    void update_cost_table(const std::vector<lingdist::AlignmentResult> &results, lingdist::CostTable &cost, CharMap &char_map)
    {
        std::unordered_map<std::string, double> row_counts;
        char_map.pair_count.clear();
        // we do not clear cost, since all its entries will be updated anyway
        // double total_count = 0.0;
        double eff_val;

        for (const auto &res : results)
        {
            for (const auto &path : res.alignments)
            {
                // This is not mentioned in Wieling(2012), but seems necessary:
                // If we find multiple optimal paths, we must make sure that pairs in each path only contribute fractionally,
                // so that this alignment contributes only 1.0 to the all counts.
                eff_val = 1.0 / static_cast<double>(path.size());
                for (const auto &char_pair : path)
                {

                    const std::string &c1 = char_pair.first;
                    const std::string &c2 = char_pair.second;
                    if (c1 != c2)
                    {
                        char_map.increment_item(c1, c2, eff_val);
                        insert_or_increment(row_counts, c1, eff_val);
                        insert_or_increment(row_counts, c2, eff_val);
                        // total_count += eff_val;
                    }
                }
            }
        }
        double max_raw_dist = -std::numeric_limits<double>::infinity(), min_raw_dist = std::numeric_limits<double>::infinity();
        std::vector<std::tuple<std::string, std::string, double>> update_queue;
        for (uint32_t char1_id = 0; char1_id < char_map.rev_idx.size(); ++char1_id)
        {
            const std::string &c1 = char_map.rev_idx.at(char1_id);
            for (uint32_t char2_id = char1_id; char2_id < char_map.rev_idx.size(); ++char2_id)
            {
                const std::string &c2 = char_map.rev_idx.at(char2_id);
                if (char_map.has(c1, c2))
                {
                    uint64_t key = char_map.s2k(c1, c2);
                    double c_xy = char_map.pair_count.at(key);
                    double c_x = row_counts.at(c1);
                    double c_y = row_counts.at(c2);

                    // $N$: total number of all character pairs aligned
                    // $c(x,y)$: count of co-occurrence of x and y
                    // $c(x),c(y)$: count of occurrence of x and y respectively

                    // $\mathrm{PMI} = \log_2 \frac{p(x,y)}{p(x)p(y)}=\log_2 \frac{c(x,y)N}{c(x)c(y)}$

                    // Thus $\mathrm{PMI} = \frac1{\ln2}[\ln \frac{c(x,y)}{c(x)c(y)}+ \ln N]$

                    // $N$ and $\ln 2$ are constants participating only in a linear transformation, which means
                    // after normalization, they do not affect the relative distances.
                    // we only need to update the cost table with - log(c(x, y) / (c(x) c(y)))

                    double raw_dist = -std::log((c_xy) / (c_x * c_y + 1e-10) + 1e-10); // /lingdist::LOG2; // log base 2
                    update_queue.emplace_back(c1, c2, raw_dist);
                    max_raw_dist = std::max(max_raw_dist, raw_dist);
                    min_raw_dist = std::min(min_raw_dist, raw_dist);
                }
                else
                {
                    if (c1 == c2)
                    {
                        cost.set_cost(c1, c2, 0.0);
                        cost.set_cost(c2, c1, 0.0);
                    }
                    else
                    {
                        cost.set_cost(c1, c2, 1.0);
                        cost.set_cost(c2, c1, 1.0);
                    }
                }
            }
        }

        // Now we only need to update values in update_queue to [0,1] based on min_raw_dist and max_raw_dist
        for (const auto &item : update_queue)
        {
            const std::string &c1 = std::get<0>(item);
            const std::string &c2 = std::get<1>(item);
            double raw_dist = std::get<2>(item);
            double norm_dist = (raw_dist - min_raw_dist) / (max_raw_dist - min_raw_dist);
            cost.set_cost(c1, c2, norm_dist);
            cost.set_cost(c2, c1, norm_dist);
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
        result.distance = nwords > 0 ? sum_dist / nwords : -1.0;
        return result;
    }
}

namespace lingdist
{
    List pmi_df(const DataFrame &data, const String &delim, bool squareform, bool parallel, int n_threads, int max_epochs, double tol, int alignment_max_paths, bool quiet)
    {
        bool verbose = !quiet;
        if (verbose)
            Rprintf("Starting PMI distance computation on data frame with %d rows...\n", data.nrow());

        // First, prepare a default cost table
        lingdist::StrVec unique_chars = lingdist::get_all_unique_syms(data, delim, true);

        if (verbose)
            Rprintf("Identified %d unique characters in the data frame.\n", static_cast<int>(unique_chars.size()));
        lingdist::CostTable prev_cost, cost = lingdist::build_default_cost_table(unique_chars);
        if (verbose)
            Rprintf("Initialized default cost table.\n");

        List report;
        CharMap char_map(unique_chars);
        auto rows_vector = lingdist::split_df(data, delim);
        auto [lab1Col, lab2Col, row_pairs] = lingdist::get_row_pairs(data, true);

        int n_row_pairs = static_cast<int>(row_pairs.size());
        if (verbose)
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
            if (verbose)
                Rprintf("Epoch %d/%d\n", iepoch, max_epochs);

            RcppThread::ProgressBar *bar = nullptr;
            if (verbose)
                bar = new RcppThread::ProgressBar(n_row_pairs, 1);
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

            if (verbose)
                Rprintf("Updating cost table...\n");
            prev_cost = cost;
            update_cost_table(results, cost, char_map);
            // Check convergence
            sum_diff = 0.0;
            for (size_t idx = 0; idx < cost.data.size(); ++idx)
            {
                sum_diff += std::fabs(cost.data[idx] - prev_cost.data[idx]);
            }
            mean_diff = sum_diff / (static_cast<double>(cost.data.size()) - static_cast<double>(unique_chars.size())); // exclude diagonal entries
            if (verbose)
                Rprintf("Cost table updated, mean absolute difference: %.6f, total absolute difference: %.6f\n", mean_diff, sum_diff);
            if (mean_diff < tol)
            {
                if (verbose)
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
        if (verbose)
            Rprintf("PMI distance computation completed.\n");
        return report;
    }
}