#include "../../include/gene_express.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <istream>
#include <numeric>
#include <spdlog/spdlog.h>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <iostream>

using namespace std;

GeneExpress::GeneExpress(string file_path) {
    read_gene_express_origin(file_path);


    string essential_protein_file = ESSENTIAL_PROTEIN;
    read_essential_protein(essential_protein_file);
    string protein_location_file = "/home/jh/code/JHPC/dataset/Yeast/DAG/subcellular1.txt";
    read_protein_location(protein_location_file);
}

GeneExpress::~GeneExpress() { gene_express.clear(); }

// 读取蛋白质的亚细胞定位
void GeneExpress::read_protein_location(string& file_path) {
    fstream file(file_path);
    if (!file.is_open()) {
        spdlog::error("Failed to open file!, {}", file_path);
        exit(1);
    }
    string line;
    string protein_name;
    string location;

    while (getline(file, line)) {
        istringstream iss(line);

        iss >> protein_name >> location;

        if (protein_location.find(location) == protein_location.end()) {
            vector<string> temp{location};
            protein_location.insert(
                std::move(make_pair(std::move(protein_name), std::move(temp))));
        } else {
            protein_location[protein_name].push_back(location);
        }
    }
    // std::cout << "<INFO> Read " << protein_location.size() << endl;
    // for(auto& it: protein_location) {
    //     std::cout << it.first << " " << it.second.size() << "\t";
    //     for(auto& it1: it.second) {
    //         std::cout << it1 << " ";
    //     }
    //     std::cout << endl;
    // }
    file.close();
}

void GeneExpress::read_gene_express(string& file_path) {
    fstream file(file_path);
    if (!file.is_open()) {
        spdlog::error("Failed to open file!, {}", file_path);
        exit(1);
    }
    string line;
    while (getline(file, line)) {
        istringstream iss(line);
        string protein_name;
        iss >> protein_name;
        vector<double> temp_express;
        double temp;
        while(iss >> temp) {
            temp_express.push_back(temp);
        }
        gene_express.insert(
            std::move(make_pair(std::move(protein_name), std::move(temp_express))));
    }

    file.close();
}

void GeneExpress::read_gene_express_origin(string& file_path) {
    fstream file(file_path);
    if (!file.is_open()) {
        spdlog::error("Failed to open file!, {}", file_path);
        exit(1);
    }
    string line;
    getline(file, line);
    getline(file, line);
    while(getline(file, line)) {
        istringstream iss(line);
        string protein_name;
        iss >> protein_name;
        iss >> protein_name;

        vector<double> temp_express;
        double express;
        iss >> express;
        while(iss >> express) {
            temp_express.push_back(express);
        }
        gene_express.insert(
            std::move(make_pair(std::move(protein_name), std::move(temp_express))));
    }
    file.close();
}

void GeneExpress::read_essential_protein(string& file_path) {
    fstream file(file_path);
    if (!file.is_open()) {
        spdlog::error("Failed to open the essential protein file!, {}", file_path);
        exit(1);
    }

    string line;
    while (getline(file, line)) {
        istringstream iss(line);
        string protein;
        iss >> protein;
        essential_proteins.insert(protein);
    }
    file.close();
}

map<string, double> GeneExpress::activate_by_three_sigma() {
    map<string, double> result;
    // 计算mean and 。。。
    for (auto& it : gene_express) {
        // update average gene expression
        double temp_sum = 0.0;
        for (auto& express : it.second) {
            temp_sum += express;
        }
        double mean = temp_sum / it.second.size();

        // update gene expression variance
        double temp_variance = 0.0;
        for (auto& express : it.second) {
            temp_variance += pow(express - mean, 2);
        }
        double variance = temp_variance / (it.second.size() - 1);
        // cout << "variance: " << variance << endl;
        // update activate threshold
        double active_threshold = active_three_sigma(mean, variance);
        result.insert(make_pair(it.first, active_threshold));
    }
    return result;
}

double GeneExpress::active_three_sigma(double mean, double varience) {
    const double f = 1.0 / (1.0 + varience);
    const double s_varience = pow(varience, 0.5);
    return mean + 3 * s_varience * (1.0 - f);
}

// need todo()!
map<string, double> GeneExpress::active_by_top() {
    map<string, double> result;
    for (auto& it : gene_express) {
        vector<double> temp_expression(it.second.begin(), it.second.end());
        sort(temp_expression.begin(), temp_expression.end());
        // return temp_expression[10];
        result.insert(make_pair(it.first, temp_expression[12]));
    }
    return result;
}

pair<double, double> GeneExpress::calculate_mean_varience(vector<double>& arr) {
    double sum = accumulate(arr.begin(), arr.end(), 0.0);
    const double mean = sum / static_cast<double>(arr.size());
    double sum_varience = 0.0;
    for (auto value : arr) {
        sum_varience += pow(mean - value, 2);
    }
    const double varience = sum_varience / static_cast<double>(arr.size());
    return make_pair(mean, varience);
}

vector<vector<bool>> GeneExpress::update_same_location(const Graph* g) {
    vector<vector<bool>> same_location(g->node_count, vector<bool>(g->node_count, false));
    for(auto& e: g->edges) {
        auto it1 = protein_location.find(g->id_protein.find(e.smllar)->second);
        auto it2 = protein_location.find(g->id_protein.find(e.bigger)->second);
        if(it1 == protein_location.end() || it2 == protein_location.end()) {
            same_location[e.smllar][e.bigger] = false;
            same_location[e.bigger][e.smllar] = false;
            continue;
        }
        set<string> common;
        set_intersection(it1->second.begin(), it1->second.end(), it2->second.begin(), it2->second.end(), inserter(common, common.begin()));
        if(!common.empty()) {
            same_location[e.smllar][e.bigger] = true;
            same_location[e.bigger][e.smllar] = true;
        } else {
            same_location[e.smllar][e.bigger] = false;
            same_location[e.bigger][e.smllar] = false;
        }
    }
    return same_location;
}

vector<Graph> GeneExpress::build_KPIN(const Graph* g) {
    // 亚细胞定位数据加入后用处不大
    vector<vector<bool>> same_location = update_same_location(g);

    int count = gene_express.begin()->second.size();
    map<string, double> active;
    for (auto p = 0; p < g->node_count; ++p) {
        // 查找基因表达
        double active_temp = 0.0;
        auto it = gene_express.find(g->id_protein.find(p)->second);
        if (it != gene_express.end()) {
            auto mean_varience = calculate_mean_varience(it->second);
            if (essential_proteins.count(g->id_protein.find(p)->second)) {
                // 是关键蛋白质，基因表达活性阈值更低一些，以保证能在生命周期中
                // 的大多数时间内保持活性
                active_temp =
                    mean_varience.first +
                    pow(mean_varience.second, 0.5) *
                        (mean_varience.second / (1.0 + mean_varience.second));
                // active_temp = mean_varience.first;
            } else { // 非关键蛋白质，具有更高的基因表达活性表达阈值
                active_temp =
                    mean_varience.first - pow(mean_varience.second, 0.5);
            }
        } else {
            active_temp = 99999;
        }
        active.insert(make_pair(g->id_protein.find(p)->second, active_temp));
        auto it1 = gene_express.find(g->id_protein.find(p)->second);
    }

    // 构建动态网络
    vector<Graph> dpins(count);
    vector<set<string>> set_proteins(count);
    vector<vector<string>> list_edges(count);
    vector<vector<double>> edge_weight(count);
    for (int i = 0; i < g->node_count; ++i) {
        // 是够含有当前蛋白质的基因表达数据？
        string protein = g->id_protein.find(i)->second;
        auto it = gene_express.find(protein);
        if (it == gene_express.end()) {
            for (int i = 0; i < count; ++i) {
                set_proteins[i].insert(protein);
            }
            continue;
        }
        for (int i = 0; i < count; ++i) {
            if (it->second[i] > active[protein] && it->second[i] > 0) {
                set_proteins[i].insert(protein);
            }
        }
    }

    // update edges
    for (auto& e : g->edges) {
        string smaller = g->id_protein.find(e.smllar)->second;
        string bigger = g->id_protein.find(e.bigger)->second;

        for (int i = 0; i < count; ++i) {
            if (set_proteins[i].count(g->id_protein.find(e.smllar)->second) &&
                set_proteins[i].count(g->id_protein.find(e.bigger)->second)) {
            // if (set_proteins[i].count(g->id_protein.find(e.smllar)->second) &&
            //     set_proteins[i].count(g->id_protein.find(e.bigger)->second) &&
            //     same_location[e.smllar][e.bigger]) {
                // std::cout << "yes\t" << g->id_protein.find(e.smllar)->second
                // << "\t" << g->id_protein.find(e.bigger)->second << endl;

                list_edges[i].emplace_back(
                    g->id_protein.find(e.smllar)->second);
                list_edges[i].emplace_back(
                    g->id_protein.find(e.bigger)->second);
                edge_weight[i].emplace_back(e.weight);
            }
        }
    }
    for (int i = 0; i < count; ++i) {
        Graph temp_graph(list_edges[i], edge_weight[i]);
        dpins[i] = temp_graph;
    }
    return dpins;
}


// vector<Graph> GeneExpress::build_KPIN_no_dag(const Graph* g) {
//     // 亚细胞定位数据加入后用处不大
//     vector<vector<bool>> same_location = update_same_location(g);

//     int count = gene_express.begin()->second.size();
//     map<string, double> active;
//     for (auto p = 0; p < g->node_count; ++p) {
//         // 查找基因表达
//         double active_temp = 0.0;
//         auto it = gene_express.find(g->id_protein.find(p)->second);
//         if (it != gene_express.end()) {
//             auto mean_varience = calculate_mean_varience(it->second);
//             if (essential_proteins.count(g->id_protein.find(p)->second)) {
//                 // 是关键蛋白质，基因表达活性阈值更低一些，以保证能在生命周期中
//                 // 的大多数时间内保持活性
//                 active_temp =
//                     mean_varience.first +
//                     pow(mean_varience.second, 0.5) *
//                         (mean_varience.second / (1.0 + mean_varience.second));
//                 // active_temp = mean_varience.first;
//             } else { // 非关键蛋白质，具有更高的基因表达活性表达阈值
//                 active_temp =
//                     mean_varience.first - pow(mean_varience.second, 0.5);
//             }
//         } else {
//             active_temp = 99999;
//         }
//         active.insert(make_pair(g->id_protein.find(p)->second, active_temp));
//         auto it1 = gene_express.find(g->id_protein.find(p)->second);
//     }

//     // 构建动态网络
//     vector<Graph> dpins(count);
//     vector<set<string>> set_proteins(count);
//     vector<vector<string>> list_edges(count);
//     vector<vector<double>> edge_weight(count);
//     for (int i = 0; i < g->node_count; ++i) {
//         // 是够含有当前蛋白质的基因表达数据？
//         string protein = g->id_protein.find(i)->second;
//         auto it = gene_express.find(protein);
//         if (it == gene_express.end()) {
//             for (int i = 0; i < count; ++i) {
//                 set_proteins[i].insert(protein);
//             }
//             continue;
//         }
//         for (int i = 0; i < count; ++i) {
//             if (it->second[i] > active[protein] && it->second[i] > 0) {
//                 set_proteins[i].insert(protein);
//             }
//         }
//     }

//     // update edges
//     for (auto& e : g->edges) {
//         string smaller = g->id_protein.find(e.smllar)->second;
//         string bigger = g->id_protein.find(e.bigger)->second;

//         for (int i = 0; i < count; ++i) {
//             if (set_proteins[i].count(g->id_protein.find(e.smllar)->second) &&
//                 set_proteins[i].count(g->id_protein.find(e.bigger)->second)) {
//             // if (set_proteins[i].count(g->id_protein.find(e.smllar)->second) &&
//             //     set_proteins[i].count(g->id_protein.find(e.bigger)->second) &&
//             //     same_location[e.smllar][e.bigger]) {
//                 // std::cout << "yes\t" << g->id_protein.find(e.smllar)->second
//                 // << "\t" << g->id_protein.find(e.bigger)->second << endl;

//                 list_edges[i].emplace_back(
//                     g->id_protein.find(e.smllar)->second);
//                 list_edges[i].emplace_back(
//                     g->id_protein.find(e.bigger)->second);
//                 edge_weight[i].emplace_back(e.weight);
//             }
//         }
//     }

//     for (int i = 0; i < count; ++i) {
//         Graph temp_graph(list_edges[i]);
//         dpins[i] = temp_graph;
//     }
//     return dpins;
// }



// vector<Graph> GeneExpress::build_KPIN1(const Graph* g) {
//     int count = 12;
//     map<string, double> active;
//     for (auto p = 0; p < g->node_count; ++p) {
//         // 查找基因表达
//         double active_temp = 0.0;
//         auto it = gene_express.find(g->id_protein.find(p)->second);
//         if (it != gene_express.end()) {
//             auto mean_varience = calculate_mean_varience(it->second);
//             if (essential_proteins.count(g->id_protein.find(p)->second)) {
//                 // 是关键蛋白质，基因表达活性阈值更低一些，以保证能在生命周期中
//                 // 的大多数时间内保持活性
//                 active_temp =
//                     mean_varience.first -
//                     pow(mean_varience.second, 0.5) *
//                         (mean_varience.second / (1.0 + mean_varience.second));
//                 // active_temp = mean_varience.first;
//             } else { // 非关键蛋白质，具有更高的基因表达活性表达阈值
//                 active_temp =
//                     mean_varience.first + pow(mean_varience.second, 0.5);
//             }
//         } else {
//             active_temp = 0.0;
//         }
//         active.insert(make_pair(g->id_protein.find(p)->second, active_temp));
//         auto it1 = gene_express.find(g->id_protein.find(p)->second);
//     }

//     // 构建动态网络
//     vector<Graph> dpins(count);
//     vector<set<string>> set_proteins(count);
//     vector<vector<string>> list_edges(count);
//     vector<vector<double>> edge_weight(count);
//     for (int i = 0; i < g->node_count; ++i) {
//         // 是够含有当前蛋白质的基因表达数据？
//         string protein = g->id_protein.find(i)->second;
//         auto it = gene_express.find(protein);
//         if (it == gene_express.end()) {
//             for (int i = 0; i < count; ++i) {
//                 set_proteins[i].insert(protein);
//             }
//             continue;
//         }
//         for (int i = 0; i < count; ++i) {
//             double exp =
//                 (it->second[i] + it->second[i + 12] + it->second[i + 24]) / 3;
//             if (exp > active[protein]) {
//                 set_proteins[i].insert(protein);
//             }
//         }
//     }

//     // update edges
//     for (auto& e : g->edges) {
//         string smaller = g->id_protein.find(e.smllar)->second;
//         string bigger = g->id_protein.find(e.bigger)->second;
//         // std::cout << smaller << "\t" << bigger << endl;
//         for (int i = 0; i < count; ++i) {
//             if (set_proteins[i].count(g->id_protein.find(e.smllar)->second) &&
//                 set_proteins[i].count(g->id_protein.find(e.bigger)->second)) {
//                 // std::cout << "yes\t" << g->id_protein.find(e.smllar)->second
//                 // << "\t" << g->id_protein.find(e.bigger)->second << endl;

//                 list_edges[i].emplace_back(
//                     g->id_protein.find(e.smllar)->second);
//                 list_edges[i].emplace_back(
//                     g->id_protein.find(e.bigger)->second);
//                 edge_weight[i].emplace_back(e.weight);
//             }
//         }
//     }

//     for (int i = 0; i < count; ++i) {
//         Graph temp_graph(list_edges[i]);
//         dpins[i] = temp_graph;
//     }
//     return dpins;
// }