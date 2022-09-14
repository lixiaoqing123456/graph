#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <bitset>
#include <queue>
#include "solution.h"

using namespace std;

#define MOD (100000007)
#define MAX_SOCRE (-1)

unordered_map<int, int> node_degree;

struct Machine {
    int id;
    int node_cost;
    int edge_cost;
    int mem;
    int comm_cost;
    int acc_cal_time;
    int acc_com_time;
};


bool cmp_node(const int& v1, const int& v2)
{
    return node_degree[v1] < node_degree[v2];
}


inline void get_total_time(
    int edge_nums,
    unordered_set<int>& batch_nodes,
    Machine& m,
    vector<Machine>& machines,
    vector<bitset<128>>& machine_sets,
    unordered_map<int, bitset<128>>& node_dist,
    /*以下为输出*/
    unsigned int& score,
    unsigned int& cal_time,
    vector<int>& comm_time,
    int& need_mem
)
{
    need_mem = 2 * edge_nums;
    int not_exist_num = 0;
    int machine_num = machine_sets.size();
    int comm_cost;

    // 计算机器存在的节点，和不存在的节点
    // 对于内存，不存在的节点需要内存
    for (auto v : batch_nodes) {
        auto& n_dist = node_dist[v];
        if ((machine_sets[m.id] & n_dist) == 0) {
            not_exist_num++;
            for (int i = 0; i < machine_num; ++i) {
                if ((machine_sets[i] & n_dist) == 0) continue;
                comm_cost = m.comm_cost + machines[i].comm_cost;
                comm_time[m.id] += comm_cost;
                comm_time[i] += comm_cost;
            }
        }
    }

    need_mem += not_exist_num;
    if (m.mem < need_mem) {
        score = MAX_SOCRE;
        return;
    }
    // 计算时间
    cal_time = not_exist_num * m.node_cost + edge_nums * m.edge_cost;

    score = cal_time + *max(comm_time.begin(), comm_time.end());
    // 对于计算，edge_cost* edge_num + not_exist * node_cost
    // 对于通信，原本存在的节点不增加通信时间，只有不存在的节点才会增加通信时间

    // 输出 机器i需要减少的内存，机器I需要增加的计算时间，机器I需要增加的通信时间，其他机器需要更新的通信时间
    // 避免重复计算
}
void insert_edges(vector<pair<int, int>*>& batch_edges,
    int len,
    unordered_set<int> batch_nodes,
    vector<Machine>& machines,
    vector<bitset<128>>& machine_sets,
    unordered_map<int, bitset<128>>& node_dist)
{
    unsigned int min_socre = MAX_SOCRE;
    int machine_id = 0;
    int machine_num = machines.size();
    for (int i = 0; i < machine_num; ++i) {

    }
}
void Solution::func(const string& input_file, const string& result_file)
{
    ifstream ifs(input_file);
    int n, m, src, dst, machine_num, machine_id;
    ifs >> n >> m;
    node_degree.reserve(n);
    vector<pair<int, int>> graph;
    graph.reserve(m);
    unordered_map<int, vector<pair<int, int>*>> adj_list;
    adj_list.reserve(n);
    vector<int> nodes;
    nodes.reserve(n);
    unordered_map<int, bitset<128>> node_dist;
    // 读取图
    for (int i = 0; i < m; ++i) {
        ifs >> src >> dst;
        graph.emplace_back(src, dst);
        if (node_degree[src] == 0) {
            nodes.emplace_back(src);
        }
        if (node_degree[dst] == 0) {
            nodes.emplace_back(dst);
        }
        node_degree[src]++;
        node_degree[dst]++;
    }
    // 读取机器
    ifs >> machine_num;
    vector<Machine> machines(machine_num);
    vector<vector<pair<int, int>*>> result;
    vector<bitset<128>> machine_sets(machine_num);
    for (int i = 0; i < machine_num; ++i) {
        ifs >> machines[i].id
            >> machines[i].node_cost
            >> machines[i].edge_cost
            >> machines[i].mem
            >> machines[i].comm_cost;
    }
    // 初始化machine set
    for (int i = 0; i < machine_num; ++i) {
        machine_sets[i].set(i);
    }
    // 构建邻接表
    for (int i = 0; i < m; ++i) {
        if (node_degree[graph[i].first] < node_degree[graph[i].second]) {
            adj_list[graph[i].first].push_back(&graph[i]);
        }
        else {
            adj_list[graph[i].second].push_back(&graph[i]);
        }
    }
    sort(nodes.begin(), nodes.end(), cmp_node);
    // 正式开始插入

    int batch_put_num = 10000;
    vector<pair<int, int>*> batch_edges(batch_put_num);
    unordered_set<int> batch_nodes;
    unordered_map<int, int> visit;
    int len;

    int que_size, node_id;
    queue<int> que;
    for (auto v : nodes) {
        if (visit[v] == 1) continue;
        que.push(v);
        while (!que.empty()) {
            que_size = que.size();
            for (int s = 0; s < que_size; ++s) {
                node_id = que.front();
                que.pop();
                visit[node_id] = 1;
                auto edges = adj_list[node_id];
                for (auto e : edges) {
                    if (len >= batch_put_num) {
                        // 选择机器并插入
                        // 清空
                        batch_edges.clear();
                        batch_nodes.clear();
                        len = 0;
                    }
                    batch_edges[len] = e;
                    ++len;
                    batch_nodes.insert(e->first);
                    batch_nodes.insert(e->second);
                    if (visit[e->first] == 0) {
                        que.push(e->first);
                    }
                    if (visit[e->second] == 0) {
                        que.push(e->second);
                    }
                }

            }
        }
    }

}