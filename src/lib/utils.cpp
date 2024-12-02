#include "../../include/utils.h"
#include <algorithm>
#include <cstddef>
#include <iterator>
#include <set>
#include <spdlog/spdlog.h>
#include <string>
#include <vector>

using namespace std;

// 从小到大排序
template <typename T> bool compare(const set<T>& a, const set<T>& b) {
    return a.size() < b.size();
}

std::vector<std::set<std::string>>
function1(std::set<std::set<std::string>> complex) {
    vector<set<string>> result;
    vector<set<string>> complex_vec(complex.begin(), complex.end());
    std::sort(complex_vec.begin(), complex_vec.end(), compare<string>);

    for (size_t i = 0; i < complex_vec.size(); ++i) {
        bool flag = true;
        for (size_t j = i + 1; j < complex_vec.size(); ++j) {
            set<string> common;
            set_intersection(complex_vec[i].begin(), complex_vec[i].end(),
                             complex_vec[j].begin(), complex_vec[j].end(),
                             std::inserter(common, common.begin()));
            if (common.size() == complex_vec[i].size()) {
                flag = false;
                break;
            }
        }
        if (flag) {
            result.push_back(complex_vec[i]);
        }
    }

    return result;
}

std::vector<std::set<std::string>> function2(std::set<std::set<std::string>> complex) {
    std::vector<std::set<std::string>> result;
    std::vector<std::set<std::string>> complex_vec(complex.begin(), complex.end());
    std::sort(complex_vec.begin(), complex_vec.end(), compare<string>);

    for (size_t i = 0; i < complex_vec.size(); ++i) {
        bool flag = true;
        for (size_t j = i + 1; j < complex_vec.size(); ++j) {
            std::set<std::string> common;
            std::set_intersection(complex_vec[i].begin(), complex_vec[i].end(),
                                  complex_vec[j].begin(), complex_vec[j].end(),
                                  std::inserter(common, common.begin()));
            size_t setSize = std::max(complex_vec[i].size(), complex_vec[j].size());
            double os_score = setSize > 0 ? static_cast<double>(common.size()) / setSize : 0.0;
            if (os_score > 0.85) {
                flag = false;
                break; // 一旦发现重叠超过阈值，跳出内部循环
            }
        }
        if (flag) {
            result.push_back(complex_vec[i]); // 如果没有超过阈值，添加到结果
        }
    }

    return result;
}


bool hasOverlap(const std::set<std::string>& complex1, const std::set<std::string>& complex2) {
    for (const std::string& protein : complex1) {
        if (complex2.find(protein) != complex2.end()) {
            // 如果发现共享蛋白质，则存在重叠
            return true;
        }
    }
    return false;
}

// 合并蛋白质复合物
std::set<std::string> mergeComplexes(const std::set<std::string>& complex1, const std::set<std::string>& complex2) {
    std::set<std::string> mergedComplex = complex1;
    for (const std::string& protein : complex2) {
        mergedComplex.insert(protein);
    }
    return mergedComplex;
}

// 删除冗余的蛋白质复合物
void deleteHonoraryComplexes(std::set<std::set<std::string>>& complexes) {
    std::vector<std::set<std::string>> toRemove; // 用来存储需要删除的复合物

    for (const auto& complex1 : complexes) {
        bool isHonorary = true; // 假设当前复合物是荣誉的
        for (const auto& complex2 : complexes) {
            if (complex1 != complex2 && hasOverlap(complex1, complex2)) {
                // 找到了重叠，所以当前复合物不是荣誉的
                isHonorary = false;
                break;
            }
        }
        if (isHonorary) {
            toRemove.push_back(complex1); // 如果复合物是荣誉的，则添加到删除列表
        }
    }

    // 删除荣誉复合物
    for (const auto& complex : toRemove) {
        complexes.erase(complex);
    }
}

