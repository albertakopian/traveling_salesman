#include<iostream>
#include<vector>
#include<map>
#include<random>
#include<algorithm>
#include<cmath>
#include<unordered_map>
#include<random>
#include<chrono>

#include "POINT.h"
#include "EDGE.h"


void point_generation(int n, std::unordered_map<int, point> &hash_points)
{
    const double PI = 3.1415926535;
    std::uniform_real_distribution<double> uniDist(0, 2 * PI);
    std::normal_distribution<double> normDist;
    std::vector<double> results;
    std::default_random_engine point_generation((unsigned)std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()));
    for (int i = 0; i < n; ++i)
    {
        double r = normDist(point_generation);
        double phi = uniDist(point_generation);
        double x = r * (cos(phi));
        double y = r * (sin(phi));
        hash_points.emplace(i, point(x, y));
    }
}

void find_mst(const int &n, std::unordered_map<int, point> &hash_points, std::vector<std::vector<int>> &graph_mst)
{
    std::vector<int> prev(n, 0);
    std::vector<int> visited(n, 0);
    std::vector<int> mst;
    double length = 0;
    std::multimap<double, std::pair<int, int>> Q;
    Q.emplace(0, std::make_pair( 0, -1 ));

    while (!Q.empty())
    {
        auto cur = Q.begin();
        Edge cur_edge = Edge(cur->second.first, cur->second.second, cur->first);
        int cur_vertex = cur_edge.first_vertex;
        if (visited[cur_vertex])
        {
            Q.erase(cur);
            continue;
        }
        visited[cur_vertex] = 1;
        prev[cur_edge.second_vertex] = cur_vertex;
        length += cur_edge.len;

        for (int i = 0; i < n; ++i)
            if (!visited[i])
                Q.emplace(distance(hash_points[cur_vertex], hash_points[i]), std::make_pair(i, cur_vertex));

        mst.push_back(cur_vertex);
        Q.erase(cur);
    }
    for (int i = 1; i < mst.size(); ++i)
    {
        graph_mst[i].push_back(prev[i]);
        graph_mst[prev[i]].push_back(i);
    }
}

void dfs(int vertex, std::vector<int> &visited, std::vector<std::vector<int>> &graph_mst, std::vector<int>& dfs_in_order){
    visited[vertex] = true;
    dfs_in_order.push_back(vertex);
    for (int next_vertex: graph_mst[vertex]){
        if (!visited[next_vertex]){
            dfs(next_vertex, visited, graph_mst, dfs_in_order);
        }
    }
}

double variant_A_answer(const int& n, std::vector<std::vector<int>> &graph_mst, std::unordered_map<int, point> &hash_points)
{
    double value = 0;
    std::vector<int> visited(n, false);
    std::vector<int> dfs_in_order;
    dfs(0, visited, graph_mst, dfs_in_order);
    dfs_in_order.push_back(0);
    for (int i = 0; i < n; ++i){
        value += distance(hash_points[dfs_in_order[i]], hash_points[dfs_in_order[i + 1]]);
    }
    return value;
}

double salesman_path(std::vector<int> &path, std::unordered_map<int, point> &hash_points)
{
    double value = 0;
    for (int i = 0; i < path.size(); ++i)
        value += distance(hash_points[path[i]], hash_points[path[(i + 1) % path.size()]]);
    return value;
}

double really_answer(const int& n, std::unordered_map<int, point>& hash_points)
{
    std::vector<int> path(n);
    for (int i = 0; i < n; ++i)
        path[i] = i;
    double ans = 1000000000.0;
    do {
        ans = std::min(ans, salesman_path(path, hash_points));
    } while (std::next_permutation(path.begin(), path.end()));
    return ans;
}

double standard_deviation(std::vector<double> &v_real, std::vector<double> &v_our, const int& num_steps)
{
    double sum_value = 0;
    for (int i = 0; i < num_steps; ++i)
        sum_value += v_our[i] / v_real [i];
    sum_value /= num_steps;
    double avg_square = 0;
    for (int i = 0; i < num_steps; ++i)
        avg_square += (v_our[i] / v_real[i] - sum_value) * (v_our[i] / v_real[i] - sum_value);
    return sqrt(avg_square / num_steps);
}

void colculate_the_result(const int&n, const int& num_steps, std::vector<std::vector<double>> &real_ans, std::vector<std::vector<double>> &variant_A_ans)
{
    for (int i = 2; i < n; ++i){
        for (int j = 0; j < num_steps; ++j){
            std::unordered_map<int, point> hash_points;
            std::vector<std::vector<int>> graph_mst(i);
            point_generation(i, hash_points);
            double r_a = really_answer(i, hash_points);

            find_mst(i, hash_points, graph_mst);
            double t_a = variant_A_answer(i, graph_mst, hash_points);
            real_ans[i].push_back(r_a);
            variant_A_ans[i].push_back(t_a);
        }
    }
}

int main()
{
    const int num_steps = 30;
    const int n = 10;

    std::vector<std::vector<double>> real_ans(n);
    std::vector<std::vector<double>> variant_A_ans(n);

    colculate_the_result(n, num_steps, real_ans, variant_A_ans);

    for (int i = 2; i < n; ++i){
        std::cout  << i << " : " << standard_deviation(real_ans[i], variant_A_ans[i], num_steps) << '\n';
    }

}