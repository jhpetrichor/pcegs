#ifndef __P_VALUE_H__
#define __P_VALUE_H__

#include <string>
#include <cmath>
#include <set>
#include <unordered_map>
#include <vector>


#define PROTEIN_GO "/home/jh/code/JHPC/dataset/Yeast/DAG/go-slim.txt"

using namespace std;

class PValue {
public:
    vector<set<string>> complexes;
    PValue(string& ppi_file, string& complex_file);
    long double calculate_p_value(set<string>& proteins) const;

private:
    unordered_map<string, set<string>> go_protein;
    unordered_map<string, set<string>> protein_go;
    set<string> proteins;

    void initialize();
    // double calculate_p_value() const;
    static unsigned long long combination(int n, int k);
    void read_ppi_porotein(string& ppi);
    void read_complex(string& complex_file);
};

#endif