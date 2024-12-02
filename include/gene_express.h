#ifndef __GENE_EXPRESS_H__
#define __GENE_EXPRESS_H__

#include "config.h"
#include "graph.h"

#include <map>
#include <utility>
#include <vector>
#include <unordered_map>

using namespace std;

enum DPIN_MEHTOD{
    THREE_SIGMA,   // three-sigma
    TOP,           // 前。。。
    TIME,          // 时间周期
    ESSENTIAL
};

enum Level {
    HIGH,
    MEDIUM,
    LOW,
};

class GeneExpress {
public:
    map<string, vector<double>> gene_express;
    map<string, vector<string>> protein_location;
    set<string> essential_proteins;
    vector<vector<bool>> same_location;


    unordered_map<string, vector<Level>> expression_level;
    unordered_map<string, pair<double, double>> level_probability;

    // 更新两个蛋白质是否有相同的亚细胞定位，结果返回
    vector<vector<bool>> update_same_location(const Graph* g);
    explicit GeneExpress(string file_path = GENE_EXPRESSION);
    ~GeneExpress();
    void read_gene_express(string& file_path);
    void read_gene_express_origin(string& file_path);
    void read_protein_location(string& file_path);
    void read_essential_protein(string& file_path);
    map<string, double> activate_by_three_sigma();
    map<string, double> active_by_top();
    // vector<UnGraph> build_KPIN(const UnGraph* g);
    vector<Graph> build_KPIN(const Graph* g);
    vector<Graph> build_KPIN_3_sigma(const Graph* g);
    vector<Graph> build_KPIN_no_dag(const Graph* g);
    vector<Graph> build_KPIN1(const Graph* g);

private:
    static double active_three_sigma(double mean, double varience);
    pair<double, double> 
    static inline calculate_mean_varience(vector<double>& arr);
};


#endif
