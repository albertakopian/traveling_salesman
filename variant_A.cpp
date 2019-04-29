
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

using::std::cin;
using::std::cout;
using::std::vector;


void generate(int n, std::unordered_map<int, point> &hash_points)
{
    const double PI = 3.1415926535;
    std::uniform_real_distribution<double> uniDist(0, 2 * PI);
    std::normal_distribution<double> normDist;
    vector<double> results;
    std::default_random_engine generator((unsigned)std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()));
    for (int i = 0; i < n; ++i)
    {
        double r = normDist(generator);
        double phi = uniDist(generator);
        double x = r * (cos(phi));
        double y = r * (sin(phi));
        hash_points.emplace(i, point(x, y));
    }
}

void Prima(const int &n, std::unordered_map<int, point> &hash_points, vector<vector<int>> &graph_mst)
{
    vector<int> prev(n, 0);
    vector<int> visited(n, 0);
    vector<int> mst;
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


void dfs(int&step, int vertex, int& k, const int& n, double&value, vector<int> &visited, vector<vector<int>> &graph_mst, std::unordered_map<int, point> &hash_points)
{
    visited[vertex] = step;
    step += 1;
    value += distance(hash_points[vertex], hash_points[k]);
    k = vertex;

    for (int i = 0; i < graph_mst[vertex].size(); ++i)
        if (visited[graph_mst[vertex][i]] == 0)
            dfs(step, graph_mst[vertex][i], k, n, value, visited, graph_mst, hash_points);
}

double triangle_MST(const int& n, vector<vector<int>> &graph_mst, std::unordered_map<int, point> &hash_points)
{
    double value = 0;
    vector<int> visited(n, 0);
    int k = 0;
    int step = 1;
    dfs(step, 0, k, n, value, visited, graph_mst, hash_points);
    for (int i = 0; i < n; ++i)
    {
        if (visited[i] == step - 1)
        {
            value += distance(hash_points[0], hash_points[i]);
            break;
        }
    }
    return value;
}

double tree_path(vector<int> &list, std::unordered_map<int, point> &hash_points)
{
    double value = 0;
    for (int i = 0; i < list.size(); ++i)
        value += distance(hash_points[list[i]], hash_points[list[(i + 1) % list.size()]]);
    return value;
}

double really_answer(const int& n, std::unordered_map<int, point>& hash_points)
{
    vector<int> list(n);
    for (int i = 0; i < n; ++i)
        list[i] = i;
    double ans = 1000000000.0;
    do {
        ans = std::min(ans, tree_path(list, hash_points));
    } while (std::next_permutation(list.begin(), list.end()));
    return ans;
}

double avg_square_deviation_quality(vector<double> &v_real, vector<double> &v_our, const int& num_steps)
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

void fill_value_of_path(const int&n, const int& num_steps, vector<vector<double>> &real_ans, vector<vector<double>> &triangle_ans)
{
    for (int i = 2; i < n; ++i){
        for (int j = 0; j < num_steps; ++j){
            std::unordered_map<int, point> hash_points;
            vector<vector<int>> graph_mst(i);
            generate(i, hash_points);
            double r_a = really_answer(i, hash_points);

            Prima(i, hash_points, graph_mst);
            double t_a = triangle_MST(i, graph_mst, hash_points);
            real_ans[i].push_back(r_a);
            triangle_ans[i].push_back(t_a);
            cout << ".";
        }
    }
}

int main()
{
    const int num_steps = 30;
    const int n = 10;

    vector<vector<double>> real_ans(n);
    vector<vector<double>> triangle_ans(n);

    fill_value_of_path(n, num_steps, real_ans, triangle_ans);

    for (int i = 2; i < n; ++i){
        cout << "for maching " << i << " : " << avg_square_deviation_quality(real_ans[i], triangle_ans[i], num_steps) << '\n';
    }

}