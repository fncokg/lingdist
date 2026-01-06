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
        std::vector<Point> curr_path;

        // Now we start from the bottom-right corner
        int curr_x = nrow - 1, curr_y = ncol - 1;
        while (true)
        {
            if (curr_x == 0 && curr_y == 0)
            {
                curr_path.insert(curr_path.begin(), std::make_pair(curr_x, curr_y));
                result.push_back(curr_path);
                if (stack.empty())
                {
                    break;
                }
                else
                {
                    auto stack_back = stack.back();
                    stack.pop_back();
                    curr_x = stack_back.first.first;
                    curr_y = stack_back.first.second;
                    curr_path = stack_back.second;
                }
            }
            curr_path.insert(curr_path.begin(), std::make_pair(curr_x, curr_y));
            Points next_points = dist[curr_x][curr_y].second;
            curr_x = next_points[0].first;
            curr_y = next_points[0].second;
            if (next_points.size() > 1)
            {
                for (std::size_t i = 1; i < next_points.size(); i++)
                {
                    stack.push_back(std::make_pair(next_points[i], curr_path));
                    if (static_cast<int>(result.size()) >= max_paths)
                    {
                        break;
                    }
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
