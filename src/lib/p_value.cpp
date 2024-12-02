#include "../../include/p_value.h"
#include <algorithm>
#include <limits>
#include <iomanip>
#include <iostream>
#include <fstream>


PValue::PValue(string& ppi_file, string& complex_file) {
    read_ppi_porotein(ppi_file);
    read_complex(complex_file);
    initialize();

}

unsigned long long PValue::combination(int n, int k){
    if (k > n - k) {
        k = n - k;  // 优化，减少计算量
    }
    unsigned long long numerator = 1;
    unsigned long long denominator = 1;
    for (int i = 0; i < k; ++i) {
        numerator *= (n - i);
        denominator *= (i + 1);
    }
    return numerator / denominator;
}


void PValue::read_ppi_porotein(string& ppi) {
    fstream ppi_file(ppi);
    if(!ppi_file.is_open()) {
        std::cerr << "Error opening file: " << ppi << std::endl;
        exit(EXIT_FAILURE);
    }
    string ppi_line;
    while(getline(ppi_file, ppi_line)) {
        stringstream ss(ppi_line);
        string protein1, protein2;
        ss >> protein1 >> protein2;
        proteins.insert(protein1);
        proteins.insert(protein2);
    }
    ppi_file.close();
}

void PValue::initialize() {
    fstream file(PROTEIN_GO);
    if(!file.is_open()) {
        std::cerr << "Error opening file: " << PROTEIN_GO << std::endl;
        exit(EXIT_FAILURE);
    }

    string line;
    while(getline(file, line)) {
        stringstream ss(line);
        string protein, go;
        ss >> protein;

        //只讨论出现在PPI网络中的蛋白质
        if(proteins.count(protein) == 0) {
            continue;
        }

        auto it = protein_go.find(protein);
        if(it == protein_go.end()) {
            protein_go[protein] = set<string>();
        }

        while(ss >> go) {
            protein_go[protein].insert(go);
            // go_protein
            auto it_go = go_protein.find(go);
            if(it_go == go_protein.end()) {
                go_protein[go] = set<string>{protein};
            } else {
                it_go->second.insert(protein);
            }
        }
    }
    file.close();
}


void PValue::read_complex(string& complex_file){
    fstream complex_f(complex_file);
    if(!complex_f.is_open()) {
        std::cerr << "Error opening file: " << complex_file << std::endl;
        exit(EXIT_FAILURE);
    }
    string line;
    while(getline(complex_f, line)) {
        set<string> complex;
        stringstream ss(line);
        string protein;
        while(ss >> protein) {
            complex.insert(protein);
        }
        complexes.emplace_back(complex);
    }
    complex_f.close();
}

long double PValue::calculate_p_value(set<string>& complex) const {
    // 用于注释该蛋白质的所有GO terms
    set<string> gos;
    for(auto protein : complex) {
        auto it = protein_go.find(protein);
        if(it == protein_go.end()) {
            continue;
        }
        gos.insert(it->second.begin(), it->second.end());
    }
    
    int N = proteins.size();
    int C = complex.size();
    
    // 对于每一个GO term
    for(auto go : gos) {
        auto it = go_protein.find(go);
        if(it != go_protein.end()) {
            int F = it->second.size();
            vector<string> common;
            set_intersection(complex.begin(), complex.end(), it->second.begin(), it->second.end(), back_inserter(common));
            int m = common.size();

            // 计算p-value
            long double p_value = 1.0;
            std::cout << "N: " << N << ", C: " << C << ", F: " << F << ", m" << m << endl;
            for(int i = 0; i < m; ++i) {
                long double term = (long double)(combination(F, i) * combination(N-F, C-i)) / (long double)combination(F, C-i);
                p_value -= term;
                // Limit the number of decimal digits to prevent precision loss
                p_value = std::round(p_value * std::pow(10, std::numeric_limits<long double>::max_digits10 - 1)) / std::pow(10, std::numeric_limits<long double>::max_digits10 - 1);
            }
            std::cout << std::setprecision(std::numeric_limits<long double>::max_digits10) << p_value << std::endl;
        }
    }
    return 0.0; // 未实现
}