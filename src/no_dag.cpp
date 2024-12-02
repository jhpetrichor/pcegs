#include "../include/graph.h"
#include "../include/config.h"
#include "../include/gene_express.h"
#include "../include/utils.h"

#include <cstdlib>
#include <filesystem>
#include <spdlog/spdlog.h>

#include <cmath>
#include <string>
#include <vector>

using namespace std;

const filesystem::path path =
    std::filesystem::current_path().parent_path() / "dataset/Yeast/PPI";
const filesystem::path result_path = std::filesystem::current_path().parent_path() / "result";

class KPIN_CNS {
  public:
    Graph* g;
    double ct;
    double at;
    
    vector<double> node_weight;        // 节点的权重
    vector<vector<double>> attraction; // 节点之间的相互吸引力
    vector<double> node_attraction;    // 节点的影响力之和
    vector<vector<double>> uni_sim;

    KPIN_CNS(Graph* _g, double _ct, double _at) {
        g = _g;
        ct = _ct;
        at = _at;
    }
    set<set<string>> find_complex() {
        calculate_node_weight();
        calculate_attraction();
        update_uni_sim();
        set<set<string>> complexes = _m_find_complex();
        return complexes;
    }

  private:
    void calculate_node_weight() {
        node_weight.resize(g->node_count, 0.0);
        for (auto& e: g->edges) {
            node_weight[e.smllar] += e.weight;
            node_weight[e.bigger] += e.weight;
        }
    }

    void calculate_attraction() {
        attraction.resize(g->node_count, vector<double>(g->node_count, 0));
        node_attraction.resize(g->node_count, 0.0);
        for (auto& e: g->edges) {
            const double dis = 1.0 - log(1.0 / (1 + e.weight));
            const double f = node_weight[e.smllar] * node_weight[e.bigger] / pow(dis, 2);
            node_attraction[e.smllar] += f;
            node_attraction[e.bigger] += f;
            attraction[e.smllar][e.bigger] = f;
            attraction[e.bigger][e.smllar] = f;
        }
    }

    void update_uni_sim() {
        uni_sim.resize(g->node_count, vector<double>(g->node_count, 0.0));
        for(auto& e: g->edges) {
            pair<double, double> it = g->unidirectional_similarity(e.smllar, e.bigger);
            uni_sim[e.smllar][e.bigger] = it.first;
            uni_sim[e.bigger][e.smllar] = it.second;
        }
    }

    set<set<string>> _m_find_complex() {
        set<set<string>> complexes;    // core？
        for(int i = 0; i < g->node_count; ++i) {
            set<int> temp_complex{i};
            set<int> neighbor = g->get_neighbor(i);
            for(auto& j: neighbor) {
                if(uni_sim[i][j] > ct) {
                    temp_complex.insert(j);
                }
            }

            // candidate
            set<int> candidate;  // attrachment
            for(auto& n: temp_complex) {
                auto neighbor = g->get_neighbor(n);
                for(auto& nei: neighbor) {
                    if(!temp_complex.count(nei) && !candidate.count(nei)) {
                        // 计算和社区内的节点的吸引力之和
                        double sum = 0.0;
                        for(auto& c: temp_complex) {
                            sum += attraction[c][nei];
                        }
                        if(sum / node_attraction[nei]  >= at) {
                            candidate.insert(nei);
                        }
                    }
                }
            }
            temp_complex.insert(candidate.begin(), candidate.end());
            if(g->calculate_complex_cohesion_score(temp_complex) < 0.5) continue;
            if(temp_complex.size() <= 2) continue;
            set<string> complex;
            for(auto& n: temp_complex) {
                complex.insert(g->id_protein[n]);
            }

            complexes.insert(complex);
        }
        return complexes;
    }
};

DAG Graph::ancestor_child(IS_A_A2C, PART_OF_A2C);
DAG Graph::child_ancestor(IS_A_C2A, PART_OF_C2A);
int main(int argc, char const** argv) {
    string ppi_path = (path / (string(argv[1]) + ".txt")).string();
    string result = (result_path / (string(argv[2]) + ".txt")).string();

    if (!filesystem::exists(ppi_path)) {
        spdlog::error("ppi file not found. {}", ppi_path);
        return 1;
    }

    double ct = atof(argv[3]);   // core threshold
    double at = atof(argv[4]);   // uni similarity threshold

    spdlog::info("ppi_name: {}", ppi_path);
    spdlog::info("alpha: {}, beta: {}", ct, at);
    
    Graph g(ppi_path);
    
    set<set<string>>  complexes;
    GeneExpress ges;
    vector<Graph> graphs = ges.build_KPIN(&g);
    spdlog::info("edpins count: {}", graphs.size());
    spdlog::info("begin detecting protein complexes ... ...");
    for(auto& pg: graphs) {
        pg.calculate_balanced_weight();
        pg.normalize_edge_weight_min_max();
        KPIN_CNS kpin(&pg, ct, at);
        auto temp_complexes = kpin.find_complex();
        complexes.insert(temp_complexes.begin(), temp_complexes.end());
    }
    spdlog::debug("complex count: {}", complexes.size());
    // Complex::write_complex(result, complexes);

    auto res = function2(complexes);
    Complex::write_complex(result, res);
    spdlog::info("complex count: {}", res.size());
    spdlog::info("write to file: {}", result);
    spdlog::info("finished processing.");
    return 0;
}
