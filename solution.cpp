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

#include<cassert>

#include "solution.h"

using namespace std;

#define MOD (100000007)
#define MAX_SOCRE (-1)
#define MIN_SCORE (-2147483648)
#define MAX_MACHINE_NUM 64

struct Machine {
    int id;
    int node_cost;
    int edge_cost;
    int mem;
    int comm_cost;
    int acc_cal_time;
    int acc_com_time;
};

unordered_map<int, int> node_degree;
unordered_map<int, bitset<MAX_MACHINE_NUM>> node_dist;
Machine* machines;
vector<bitset<MAX_MACHINE_NUM>> machine_bits;

long* g_com_cost = nullptr;
long* g_com_time = nullptr; // 表示每个机器的通信时间
long* g_com_tmp = nullptr; // 表示每条边分配到当前节点，其他节点的通信时间
long g_max_cal_time = 0;
long g_max_com_time = 0;

// 在某个范围内选取
static inline int get_best_machine(const pair<int, int>* edge,
    const bitset<MAX_MACHINE_NUM>& search_range,
    int machine_num,
    int& machine_id,
    long& max_com_time)
{
    int src = edge->first;
    int dst = edge->second;
    auto src_dist = node_dist[src];
    auto dst_dist = node_dist[dst];
    bool src_not_exist = false;
    bool dst_not_exist = false;
    bitset<MAX_MACHINE_NUM> not_search = search_range;

    int idx = 0, j = 0, row;
    long comm_cost;
    int need_mem;

    int cur_cal_time = 0;
    long cur_max_com_time = 0;
    int cur_cal_score = 0;
    int cur_com_score = 0;

    int max_score = MIN_SCORE;
    int cur_score;
    // 初始化输出
    machine_id = 0;
    max_com_time = g_max_com_time;

    while (not_search.any() && idx < machine_num) {
        if (not_search[0] == 1) { // 说明要搜索当前的机器
            // 先计算内存, 先假设存在当前机器
            src_not_exist = false;
            dst_not_exist = false;
            need_mem = 2;

            if ((src_dist & machine_bits[idx]).none()) {
                need_mem++;
                src_not_exist = true;
            }
            if ((dst_dist & machine_bits[idx]).none()) {
                need_mem++;
                dst_not_exist = true;
            }
            if (machines[idx].mem < need_mem) {
                // 如果内存不足，开始下一次循环
                not_search >>= 1;
                idx++;
                continue; // 继续下一个
            }

            // 计算时间
            cur_cal_time = machines[idx].acc_cal_time + machines[idx].edge_cost;

            // 优化两个都存在的
            if ((!src_not_exist) && (!dst_not_exist)) {
                // 两个顶点都存在
                cur_score = g_max_cal_time - cur_cal_time;
                if (cur_score > max_score) {
                    max_score = cur_score;
                    machine_id = idx;
                    max_com_time = g_max_com_time;
                }
                not_search >>= 1;
                idx++;
                continue;// 继续下一个
            }

            row = idx * machine_num;
            // 拷贝已有的通信时间
            memcpy((char*)g_com_tmp + row * sizeof(long), (char*)g_com_time, machine_num * sizeof(long));

            // 先处理src
            if (src_not_exist) { // src不在当前机器上
                // 增加计算
                cur_cal_time += machines[idx].node_cost;

                // 计算通信时间
                j = 0;
                while (src_dist.any()) {
                    if (src_dist[0] == 1) {
                        // 增加通信时间
                        comm_cost = g_com_cost[row + j];
                        g_com_tmp[row + idx] += comm_cost;
                        g_com_tmp[row + j] += comm_cost;
                    }
                    src_dist >>= 1;
                    ++j;
                }
            }

            // 处理dst
            if (dst_not_exist) { // dst不在当前机器上
                // 增加计算
                cur_cal_time += machines[idx].node_cost;
                // 计算通信时间
                j = 0;
                while (dst_dist.any()) {
                    if (dst_dist[0] == 1) {
                        // 增加通信时间
                        comm_cost = g_com_cost[row + j];
                        g_com_tmp[row + idx] += comm_cost;
                        g_com_tmp[row + j] += comm_cost;
                    }
                    dst_dist >>= 1;
                    ++j;
                }
            }

            cur_cal_score = g_max_cal_time - cur_cal_time;
            cur_max_com_time = *max_element(g_com_tmp + row, g_com_time + row + machine_num);
            cur_com_score = g_max_com_time - cur_max_com_time;
            // 计算打散到idx机器上的分数
            cur_score = cur_cal_score + cur_com_score;
            if (cur_score > max_score) {
                max_score = cur_score;
                machine_id = idx;
                max_com_time = cur_max_com_time;
            }
        } // if (not_search[0] == 1)
        not_search >>= 1;
        idx++;
    }

    return max_score;
}


static inline void put_one_edge(const pair<int, int>* edge,
    int machine_id,
    int machine_num,
    bool src_not_exist,
    bool dst_not_exist,
    long max_com_time,
    vector<pair<int, int>*>& result)
{
    // 1.插入边
    result.emplace_back(edge);

    // 2.更新计算时间和更新内存
    machines[machine_id].acc_cal_time += machines[machine_id].edge_cost;
    machines[machine_id].mem -= 2;
    if (src_not_exist) {
        machines[machine_id].acc_cal_time += machines[machine_id].node_cost;
        machines[machine_id].mem -= 1;
    }
    if (dst_not_exist) {
        machines[machine_id].acc_cal_time += machines[machine_id].node_cost;
        machines[machine_id].mem -= 1;
    }

    // 3.更新系统的通信时间
    memcpy((char*)g_com_time, (char*)g_com_tmp + machine_id * machine_num * sizeof(long), machine_num * sizeof(long));

    // 4.更新系统的最大时间
    if (g_max_cal_time < machines[machine_id].acc_cal_time) {
        g_max_cal_time = machines[machine_id].acc_cal_time;
    }
    if (g_max_com_time < max_com_time) {
        g_max_com_time = max_com_time;
    }
}

void init_com(int machine_num)
{
    auto size = sizeof(long) * machine_num * machine_num;
    int com_cost;
    g_com_cost = (long*)malloc(size);
    g_com_time = (long*)malloc(sizeof(long) * machine_num);
    g_com_tmp = (long*)malloc(size);
    memset(g_com_cost, 0, size);
    memset(g_com_time, 0, sizeof(long) * machine_num);
    memset(g_com_tmp, 0, size);
    for (int i = 0; i < machine_num; ++i) {
        for (int j = i + 1; j < machine_num; ++j) {
            com_cost = machines[i].comm_cost + machines[j].comm_cost;
            g_com_cost[i * machine_num + j] = com_cost;
            g_com_cost[j * machine_num + i] = com_cost;
        }
    }
}

bool cmp_node(const int& v1, const int& v2)
{
    return node_degree[v1] < node_degree[v2];
}

void Solution::func(const string& input_file, const string& result_file)
{
    ifstream ifs(input_file);
    int n, m, src, dst, machine_num, machine_id;
    ifs >> n >> m;

    vector<pair<int, int>> graph;
    unordered_map<int, vector<pair<int, int>*>> adj_list;
    vector<int> nodes;
    graph.reserve(m);
    adj_list.reserve(n);
    nodes.reserve(n);
    node_degree.reserve(n);
    node_dist.reserve(n);

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
    machines = (Machine*)malloc(sizeof(Machine) * machine_num);
    memset(machines, 0, sizeof(Machine) * machine_num);
    vector<vector<pair<int, int>*>> result;
    for (int i = 0; i < machine_num; ++i) {
        ifs >> machines[i].id
            >> machines[i].node_cost
            >> machines[i].edge_cost
            >> machines[i].mem
            >> machines[i].comm_cost;
    }

    // 初始化machine set
    for (int i = 0; i < machine_num; ++i) {
        machine_bits.emplace_back(bitset<MAX_MACHINE_NUM>(1 << i));
    }
    // 初始化通信相关的东西
    init_com(machine_num);

#ifdef PRINT_STAGE1
    // 阶段1 数据的输入，以及初始化，验证
    cout << "n:" << n << " nodes.size():" << nodes.size() << endl;
    cout << "m:" << m << " graph.size():" << graph.size() << endl;
    assert(n == nodes.size());
    assert(m == graph.size());
    // 打印顶点，以及度
    cout << "nodes" << endl;
    for (auto& n : nodes) {
        cout << "node:" << n << " degree:" << node_degree[n] << endl;
    }
    // 打印机器信息
    cout << "machine info" << endl;
    for (int i = 0; i < machine_num; ++i) {
        cout << machines[i].id << " " << machines[i].node_cost << " "
            << machines[i].edge_cost << " " << machines[i].mem << " " << machines[i].comm_cost << endl;
    }
    // 打印machine bits
    cout << "machine bits" << endl;
    for (int i = 0; i < machine_num; ++i) {
        cout << "id:" << machines[i].id << " bit:" << machine_bits[i] << endl;
    }
    // 打印通信相关信息
    cout << "comm info" << endl;
    for (int i = 0; i < machine_num; ++i) {
        for (int j = 0; j < machine_num; ++j) {
            cout << g_com_cost[i * machine_num + j] << " ";
        }
        cout << endl;
    }
#endif

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

    // 生成广度优先序列
    vector<pair<int, int>*> p_edges;
    p_edges.reserve(m);
    unordered_map<int, int> visit;
    visit.reserve(n);

    int node_id;
    queue<int> que;
    for (auto& v : nodes) {
        if (visit[v] == 1) continue;
        que.push(v);
        visit[v] = 1; // 表示节点曾经在队列里面
        while (!que.empty()) {
            node_id = que.front();
            que.pop();
            const auto& edges = adj_list[node_id];
            for (const auto& e : edges) {
                p_edges.emplace_back(e);
                if (visit[e->first] == 0) {
                    que.push(e->first);
                    visit[e->first] = 1;
                }
                if (visit[e->second] == 0) {
                    que.push(e->second);
                    visit[e->second] = 1;
                }
            }
        }
    }

#ifdef PRINT_STAGE2
    cout << "m:" << m << " p_edges.size()" << p_edges.size() << endl;
    assert(m == p_edges.size());
#endif

    // 搜索结果
    bitset<MAX_MACHINE_NUM> search;
    bitset<MAX_MACHINE_NUM> b1;
    bitset<MAX_MACHINE_NUM> b2;
    int max_score = MIN_SCORE;
    long max_com_time = g_max_com_time;
    for (auto& e : p_edges) {
        src = e->first;
        dst = e->second;
        const auto& src_dist = node_dist[src];
        const auto& dst_dist = node_dist[dst];
        b1 = (src_dist & dst_dist); // 不产生通信时间
        b2 = (src_dist | dst_dist);
        max_com_time = g_max_com_time;
        machine_id = 0;
        search.reset();

        if (b1.any()) { // 先搜b1
            max_score = get_best_machine(e, b1, machine_num, machine_id, max_com_time);
        }
        if (max_score > MIN_SCORE) {
            // put 
            put_one_edge(e, machine_id, machine_num,
                (src_dist & machine_bits[machine_id]).none(),
                (dst_dist & machine_bits[machine_id]).none(),
                max_com_time, result[machine_id]);
            continue;
        }
        search |= b1;
        b2 = (b2 & (~search));
        if (b2.any()) {
            // search
            max_score = get_best_machine(e, b2, machine_num, machine_id, max_com_time);

        }
        if (max_score > MIN_SCORE) {
            // put
            put_one_edge(e, machine_id, machine_num,
                (src_dist & machine_bits[machine_id]).none(),
                (dst_dist & machine_bits[machine_id]).none(),
                max_com_time, result[machine_id]);
            continue;
        }
        search |= b2;
        max_score = get_best_machine(e, ~search, machine_num, machine_id, max_com_time);
        // put
        put_one_edge(e, machine_id, machine_num,
            (src_dist & machine_bits[machine_id]).none(),
            (dst_dist & machine_bits[machine_id]).none(),
            max_com_time, result[machine_id]);
    }
}