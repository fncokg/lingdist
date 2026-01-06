#include "alignment.hpp"
#include "cost_table.hpp"

using namespace Rcpp;

namespace lingdist
{
    // dist mat -> id pairs
    // DFS to enumerate all matched optimal paths
    std::vector<std::vector<Point>> get_matched_paths(std::vector<std::vector<DistPts>> dist, int max_paths)
    {
        int nrow = static_cast<int>(dist.size());
        int ncol = static_cast<int>(dist[0].size());

        std::vector<std::vector<Point>> result;
        // A stack to store the points to explore and a path from start to the point
        std::vector<std::pair<Point, std::vector<Point>>> stack;

        // Now we start from the bottom-right corner
        stack.push_back(std::make_pair(std::make_pair(nrow - 1, ncol - 1), std::vector<Point>()));
        while (!stack.empty())
        {
            // Now we explore the top point in the stack
            int curr_x = stack.back().first.first;
            int curr_y = stack.back().first.second;
            std::vector<Point> curr_path = stack.back().second;
            stack.pop_back();

            // Add this point to the current path
            // NOTE: we add point from the end to start, so later we need to reverse the path
            // We do not insert at the beginning for performance reason: inserting at the beginning of a vector is O(n) while push_back is O(1)
            curr_path.push_back(std::make_pair(curr_x, curr_y));

            // Check if we reached the start point
            if (curr_x == 0 && curr_y == 0)
            {
                // We found a complete path, reverse it and add to result
                std::reverse(curr_path.begin(), curr_path.end());
                result.push_back(curr_path);

                // Check if we reached max paths, if so, we stop
                if (static_cast<int>(result.size()) >= max_paths)
                {
                    break;
                }
            }
            else
            {
                // We are not at the start point, so the current point must have source points
                Points prev_points = dist[curr_x][curr_y].second;

                // Push those source points to the stack for further exploration
                // The next iteration will explore the last pushed point (DFS)
                for (auto &pt : prev_points)
                {
                    // pt is a source point leading to (curr_x, curr_y)
                    // We will explore it
                    stack.push_back(std::make_pair(pt, curr_path));
                }
            }
        }
        return result;
    }

    // id pairs -> alignment result with strings
    AlignmentResult get_string_alignment_result(const StrVec &chars1, const StrVec &chars2, const CostTable &cost, int max_paths)
    {
        AlignmentResult report;
        std::string char1, char2, opt;
        int curr_x, curr_y, prev_x, prev_y;
        auto dist = edit_dist_core_dp_record(chars1, chars2, cost); // Original line
        auto paths = get_matched_paths(dist, max_paths);
        double result = dist[chars2.size()][chars1.size()].first;
        report.distance = result;

        for (auto &path : paths)
        {
            std::vector<std::pair<std::string, std::string>> alignments;
            for (std::size_t j = 1; j < path.size(); j++)
            {
                curr_x = path[j].first;
                curr_y = path[j].second;
                prev_x = path[j - 1].first;
                prev_y = path[j - 1].second;
                char1 = curr_y == prev_y ? EMPTY : chars1[curr_y - 1];
                char2 = curr_x == prev_x ? EMPTY : chars2[curr_x - 1];
                alignments.emplace_back(std::make_pair(char1, char2));
            }
            report.alignments.emplace_back(alignments);
        }
        return report;
    }

    // format alignment result as R List
    List get_string_alignment(const StrVec &chars1, const StrVec &chars2, const CostTable &cost, int max_paths)
    {
        List report;
        List path_dfs;
        AlignmentResult result = get_string_alignment_result(chars1, chars2, cost, max_paths);
        std::string char1, char2, opt;

        for (auto &path : result.alignments)
        {
            StringVector chars1_col, chars2_col, operation_col;
            NumericVector cost_col, cumcost_col;
            for (std::size_t j = 0; j < path.size(); j++)
            {
                char1 = path[j].first;
                char2 = path[j].second;
                opt = char1 == EMPTY ? "insert" : (char2 == EMPTY ? "delete" : (char1 == char2 ? "same" : "substitute"));
                chars1_col.push_back(char1);
                chars2_col.push_back(char2);
                operation_col.push_back(opt);
                cost_col.push_back(cost.get_cost(char1, char2)); // Original line
            }
            DataFrame path_df = DataFrame::create(Named("Chars1") = chars1_col, Named("Chars2") = chars2_col, Named("Operation") = operation_col, Named("Cost") = cost_col);
            path_dfs.push_back(path_df);
        }
        report["alignments"] = path_dfs;
        return report;
    }

} // namespace lingdist
