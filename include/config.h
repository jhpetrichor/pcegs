#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <filesystem>
#include <string>

const std::filesystem::path DATA_PATH = std::filesystem::current_path().parent_path() / "dataset";
const std::filesystem::path RESULT_PATH =
    std::filesystem::current_path().parent_path() / "result";

const std::string GENE_EXPRESSION = (DATA_PATH/"Yeast/DAG/gene-expression.txt").string();
const std::string ESSENTIAL_PROTEIN = (DATA_PATH/"Yeast/essential_protein.txt").string();
const std::string IS_A_A2C = (DATA_PATH/"Yeast/DAG/is_a.txt").string();
const std::string PART_OF_A2C = (DATA_PATH/"/Yeast/DAG/part_of.txt").string();
const std::string PATH_GO_TERMS = (DATA_PATH /"Yeast/DAG/go_term.txt").string();
const std::string IS_A_C2A = (DATA_PATH/"Yeast/DAG/is_a_child_ancestor.txt").string();
const std::string PART_OF_C2A = (DATA_PATH / "Yeast/DAG/part_of_child_ancestor.txt").string();

const std::string GO_SLIM = (DATA_PATH/"Yeast/DAG/go-slim.txt").string();

#define COMPLEX_MAX_SIZE 20
#endif // __CONFIG_H__
