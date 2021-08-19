#ifndef EIGENDATA_HPP
#define EIGENDATA_HPP

#include <vector>
#include <unordered_map>
#include <set>
#include <regex>

#include "setOptions.hpp"
#include "mapID.hpp"
#include "htsWrappers.hpp"
#include "dataParser.hpp"
#include "genotypeData.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

void update_blocks( int i, int j, std::vector<int>& cluster_ids, std::vector<std::set<int>>& clusters);

void read_eigen(const std::string& file_name, Eigen::SparseMatrix<double>& eig_vec, Eigen::VectorXd& eig_val, const std::vector<std::string>& kp_ids);

void write_eigen(const std::string& file_prefix, Eigen::SparseMatrix<double>& eig_vec, Eigen::VectorXd& eig_val, const std::vector<std::string>& kp_ids);

void read_sparse_GRM(const std::string& filename, Eigen::SparseMatrix<double>& GRM, const std::vector<std::string>& kp_ids, const double& r_scale, const int& r_col, std::vector<std::vector<int>>& related);

void read_dense_GRM(const std::string&, Eigen::MatrixXd&, std::vector<std::string>&);

#endif 