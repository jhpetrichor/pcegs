#include "../../include/dag.h"
#include "../../include/config.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <spdlog/spdlog.h>
#include <sstream>
#include <stack>
#include <fstream>

using namespace std;

void Node::add(Node* _ancestor) { ancestor.insert(_ancestor);}

void Node::print() const {
    spdlog::info("{} : {}", go, ancestor.size());
    for (auto _ancestor : ancestor) {
        spdlog::info("{}\t", _ancestor->go);
    }
}

/**
 * @brief:
 * @param {string} is_a_file is_a路径
 * @param {string} part_of_file part_of路径
 * @return {*}
 */
DAG::DAG(string is_a_file, string part_of_file) {
    set<string> go_terms;
    read_all_go_terms(PATH_GO_TERMS, go_terms);
    nodes.resize(go_terms.size() + 1);
    // 初始化节点
    int idx = 1;
    for (auto& g : go_terms) {
        Node* newNode = new Node(idx, g);
        GO2ID.insert(make_pair(g, idx));
        nodes[idx] = newNode;
        idx += 1;
    }
    relation.resize(go_terms.size() + 1,
                    vector<Relation>(go_terms.size() + 1, Relation::NONE));

    // add edges
    vector<vector<string>> is_a;
    vector<vector<string>> part_of;
    read_ancestor_child(is_a_file, is_a);
    read_ancestor_child(part_of_file, part_of);
    add_edges(is_a, Relation::IS_A);
    add_edges(part_of, Relation::PART_OF);

    // read protein-gos
    read_protein_go(GO_SLIM);
}

DAG::~DAG() {
    for (auto& node : nodes) {
        delete node;
    }
}

Node* DAG::addNode(int id, string& go) {
    Node* newNode = new Node(id, go);
    nodes.emplace_back(newNode);
    return newNode;
}

// 添加一条边
void DAG::addEdge(int from, int to, Relation _relation) {
    if (from == 0 || to == 0)
        return;
    Node* fromNode = nodes[from];
    Node* toNode = nodes[to];

    if (fromNode && toNode) {
        fromNode->add(toNode);
        relation[from][to] = _relation;
        relation[to][from] = _relation;
    } else {
        spdlog::error("{} {} Node not found!", from, to);
        // cout << from << " " << to << endl;
        // cerr << "Node not found!" << endl;
    }
}

// 批量添加
void DAG::add_edges(vector<vector<string>>& ancestor_child,
                    Relation _relation) {
    for (auto& item : ancestor_child) {
        if (item.size() < 2)
            continue;
        for (int i = 1; i < item.size(); ++i) {
            if (GO2ID[item[i]] == 0) {
                continue;
            }
            addEdge(GO2ID[item[i]], GO2ID[item[0]], _relation);
        }
    }
}

void DAG::print() const {
    for (auto& node : nodes) {
        node->print();
    }
}

// 计算两个节点的相似度
// double DAG::similarity(int a, int b) {
//    if (a == b) return 1;
//
//    set<Node*> aChildren = nodes[a]->ancestor;
//    set<Node*> bChildren = nodes[b]->ancestor;
//
//    int commonChildren = 0;
//    for (Node* child : aChildren) {
//        if (bChildren.find(child) != bChildren.end()) {
//            commonChildren++;
//        }
//    }
//
//    int totalChildren = aChildren.size() + bChildren.size();
//    if (totalChildren == 0) return 0.0; // 避免除以零
//
//    return static_cast<double>(2 * commonChildren) / totalChildren;
//}

// 根据相互关系添加权重，更新直接祖先的s_value
void DAG::calculate_SValue(Node* node, map<Node*, double>& s_value) {
    if (!node) {
        return;
    }
    for (auto& _ancestor : node->ancestor) {
        if (!s_value.count(_ancestor)) {
            continue;
        }
        double temp_s_value = 0.0;
        switch (relation[node->id][_ancestor->id]) {
        case Relation::IS_A:
            temp_s_value = s_value[node] * 0.8;
            break;
        case Relation::PART_OF:
            temp_s_value = s_value[node] * 0.6;
            break;
        case Relation::NONE:
            temp_s_value = 0.0;
        }
        s_value[_ancestor] = max(temp_s_value, s_value[_ancestor]);
    }
}

void DAG::find_all_paths(Node* a, Node* b, vector<Node*>& path,
                         vector<vector<Node*>>& all_path) {
    path.emplace_back(a);
    if (a == b) {
        all_path.emplace_back(path);
    } else {
        for (Node* _ancestor : a->ancestor) {
            // 祖先节点不在当前路径中
            if (!count(path.begin(), path.end(), _ancestor)) {
                find_all_paths(_ancestor, b, path, all_path);
            }
        }
    }

    // 回溯时移除当前节点
    path.pop_back();
}

void DAG::print_path(const vector<vector<Node*>>& all_path) {
    for (auto& p : all_path) {
        cout << p[0]->go;
        for (int i = 0; i < p.size() - 1; ++i) {
            cout << " --> " << p[i]->go;
        }
        cout << " --> " << p[p.size() - 1]->go << endl;
    }
}

double DAG::similarity(int a, int b) {
    // 如果是同一GO term相似性为 1
    if (a == b)
        return 1;

    // 分别获取两个GO term的共同祖先
    Node* node_a = nodes[a];
    Node* node_b = nodes[b];
    set<Node*> ancestor_a;
    set<Node*> ancestor_b;
    get_all_ancestors(node_a, ancestor_a);
    get_all_ancestors(node_b, ancestor_b);
    set<Node*> common_ancestor;
    set_intersection(ancestor_a.begin(), ancestor_a.end(), ancestor_b.begin(),
                     ancestor_b.end(),
                     inserter(common_ancestor, common_ancestor.begin()));
    // 没有共同祖先
    // 即 a 不是 b 的祖先， b也不是a的祖先，同时, a, b没有共同祖先
    if (common_ancestor.empty() && !ancestor_a.count(node_b) &&
        !ancestor_b.count(node_a)) {
        return 0.0;
    }

    set<Node*> all_ancestor;
    set_union(ancestor_a.begin(), ancestor_a.end(), ancestor_b.begin(),
              ancestor_b.end(), inserter(all_ancestor, all_ancestor.begin()));

    // 初始化 all_S-value
    map<Node*, double> all_S_value;
    for (auto& item : all_ancestor) {
        all_S_value[item] = 0.0;
    }

    // 需要迭代更新, 需分情况讨论
    if (ancestor_b.count(node_a)) { // a 是 b的祖先
        common_ancestor.insert(node_a);
        all_S_value[node_b] = 1.0;
        calculate_SValue(node_b, all_S_value);
        stack<Node*> stack;
        stack.push(node_b);
        if (!stack.empty()) {
            Node* current = stack.top();
            stack.pop();
            calculate_SValue(current, all_S_value);
            for (auto& current_ancestor : current->ancestor) {
                stack.push(current_ancestor);
            }
        }
        // 需要查找a--->b的所有路线
        vector<Node*> path;
        vector<vector<Node*>> all_path;
        find_all_paths(node_b, node_a, path, all_path);
        // 计算公共
        double sum_1 = 0.0; // 记录分子
        for (auto& common_go : common_ancestor) {
            sum_1 += all_S_value[common_go];
        }
        double sum_2 = sum_1; // 分母
        for (auto& _path : all_path) {
            for (int i = 0; i < _path.size() - 1; ++i) {
                sum_2 += all_S_value[_path[i]];
            }
        }
        return sum_1 / sum_2;
    } else if (ancestor_a.count(node_b)) { // b 是 a的祖先
        common_ancestor.insert(node_b);
        all_S_value[node_a] = 1.0;
        calculate_SValue(node_a, all_S_value);
        stack<Node*> stack;
        stack.push(node_a);
        if (!stack.empty()) {
            Node* current = stack.top();
            stack.pop();
            calculate_SValue(current, all_S_value);
            for (auto& current_ancestor : current->ancestor) {
                stack.push(current_ancestor);
            }
        }

        // 需要查找a--->b的所有路线
        vector<Node*> path;
        vector<vector<Node*>> all_path;
        find_all_paths(node_a, node_b, path, all_path);
        // 计算公共
        double sum_1 = 0.0; // 记录分子
        for (auto& common_go : common_ancestor) {
            sum_1 += all_S_value[common_go];
        }
        double sum_2 = sum_1; // 分母
        for (auto& _path : all_path) {
            for (int i = 0; i < _path.size() - 1; ++i) {
                sum_2 += all_S_value[_path[i]];
            }
        }
        return sum_1 / sum_2;
    } else { // a 和 b不清楚什么关系（或者说DAG中，ab没有联通）
        all_S_value[node_b] = 1.0;
        all_S_value[node_a] = 1.0;
        calculate_SValue(node_a, all_S_value);
        calculate_SValue(node_b, all_S_value);
        stack<Node*> stack;
        stack.push(node_a);
        stack.push(node_b);
        if (!stack.empty()) {
            Node* current = stack.top();
            stack.pop();
            calculate_SValue(current, all_S_value);
            for (auto& current_ancestor : current->ancestor) {
                stack.push(current_ancestor);
            }
        }
        double sum_1 = 0.0;
        for (auto& common_go : common_ancestor) {
            sum_1 += all_S_value[common_go];
        }
        double sum_2 = sum_1;
        vector<Node*> path;
        vector<vector<Node*>> all_path;
        find_all_paths(node_a, node_b, path, all_path);

        set<Node*> s1; // a的祖先但不是b的祖先
        set<Node*> s2; // b的祖先但不是a的祖先
        set_difference(ancestor_a.begin(), ancestor_a.end(), ancestor_b.begin(),
                       ancestor_b.end(), inserter(s1, s1.begin()));
        set_difference(ancestor_b.begin(), ancestor_b.end(), ancestor_a.begin(),
                       ancestor_a.end(), inserter(s2, s2.begin()));
        s1.insert(node_a);
        s2.insert(node_b);
        set<Node*> ss; // 表示通往共同祖先所经历的祖先
        for (auto& ii : s1) {
            for (auto& jj : common_ancestor) {
                vector<Node*> path1;
                vector<vector<Node*>> all_paths1;
                find_all_paths(ii, jj, path1, all_paths1);
                for (auto& p : all_paths1) {
                    ss.insert(p.begin(), p.end() - 1);
                }
            }
        }

        for (auto& ii : s2) {
            for (auto& jj : common_ancestor) {
                vector<Node*> path1;
                vector<vector<Node*>> all_paths1;
                find_all_paths(ii, jj, path1, all_paths1);
                for (auto& p : all_paths1) {
                    ss.insert(p.begin(), p.end() - 1);
                }
            }
        }
        for (auto& p : ss) {
            sum_2 += all_S_value[p];
        }

        return sum_1 / sum_2;
    }
}

double DAG::similarity(const string& go1, const string& go2) {
    auto it = Similarity.find(pair<string, string>{go1, go2});
    if (it != Similarity.end()) {
        return it->second;
    } else {
        double sim = 0.0;
        if (go1 == go2) {
            sim = 0.0;
        } else {
            sim = similarity(GO2ID[go1], GO2ID[go2]);
        }
        Similarity.insert(make_pair(pair<string, string>{go1, go2}, sim));
        Similarity.insert(make_pair(pair<string, string>{go2, go1}, sim));
        // if (isnan(sim) || sim <= 0.00001) {
        //     return 0.1;
        // }
        return sim;
    }
}

void DAG::read_ancestor_child(const string& file_path,
                              vector<vector<string>>& ancestor_child) {
    fstream file(file_path);
    if (!file.is_open()) {
        cerr << "Failed to open file! " << file_path << endl;
    }
    string line;
    while (getline(file, line)) {
        vector<string> _ancestor_child;
        istringstream iss(line);
        string go_term;
        while (iss >> go_term) {
            _ancestor_child.emplace_back(go_term);
        }
        ancestor_child.emplace_back(_ancestor_child);
    }
    file.close();
}

void DAG::read_all_go_terms(const string& file_path, set<string>& go_terms) {
    fstream file(file_path);
    if (!file.is_open()) {
        cerr << "Node not found!" << endl;
    }
    string line;
    while (getline(file, line)) {
        istringstream iss(line);
        string go_term;
        while (iss >> go_term) {
            go_terms.insert(go_term);
        }
    }
    file.close();
}

void DAG::get_all_ancestors(Node* node, set<Node*>& all_ancestor) {
    stack<Node*> stack;

    if (node != nullptr) {
        stack.push(node);
        while (!stack.empty()) {
            Node* current = stack.top();
            stack.pop();
            for (Node* ancestor : current->ancestor) {
                if (all_ancestor.insert(ancestor)
                        .second) {        // Insert only if it's a new ancestor
                    stack.push(ancestor); // Push the ancestor onto the stack to
                                          // explore its ancestors
                }
            }
        }
    }
}

void write_go_term(const DAG& dag) {
    ofstream file("./go_term_go.txt");
    for (auto& it : dag.GO2ID) {
        file << it.first << "\t" << it.second << endl;
    }
}

void write_nodes(const DAG& dag) {
    ofstream file("./go_term_nodes.txt");
    for (auto& it : dag.nodes) {
        file << it->go << "\t" << it->id << endl;
    }
}

// 计算两个蛋白质的相互作用, 需要传入两个蛋白质的GO term
double DAG::get_similarity_protein(const string& protein_1,
                                   const string& protein_2) {
    auto it1 = protein2gos.find(protein_1);
    auto it2 = protein2gos.find(protein_2);
    if (it1 == protein2gos.end() || it2 == protein2gos.end()) {
        return 0;
    }
    auto gos1 = it1->second;
    auto gos2 = it2->second;
    double sum_sim = 0.0;
    for (auto& g1 : gos1) {
        sum_sim += get_similarity_go_gos_by_max(g1, gos2);
    }
    for (auto& g2 : gos2) {
        sum_sim += get_similarity_go_gos_by_max(g2, gos1);
    }
    return sum_sim / (gos1.size() + gos2.size());
}

double DAG::get_similarity_go_gos_by_max(const string& go,
                                         const set<string>& gos) {
    double sim = 0.0;
    for (auto& g : gos) {
        double temp_sim = similarity(go, g);
        if (temp_sim > sim) {
            sim = temp_sim;
        }
    }
    return sim;
}

void DAG::read_protein_go(string file_path) {
    fstream file(file_path);
    string line;
    while (getline(file, line)) {
        istringstream iss(line);
        string protein;
        iss >> protein;
        set<string> go;
        string g;
        while (iss >> g) {
            go.insert(g);
        }
        protein2gos.insert(make_pair(protein, go));
    }
    file.close();
}