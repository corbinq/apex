# 21 "yax/src/readTable.hpp.c"
#ifndef READTABLE_HPP
#define READTABLE_HPP 

#include <vector>

#include "setOptions.hpp"
#include "dataParser.hpp"
#include "mapID.hpp"

#include <Eigen/Dense>


class table
{
 public:
  int header_line;

  int id_column;
  std::string id_column_name;

  int n_rows;
  int n_cols;

  id_map rows;
  id_map cols;

  Eigen::MatrixXd data_matrix;

  void setRows(std::vector<std::string>&);
  void setCols(std::vector<std::string>&);
  void readFile(const char*);
  void readHeader(const char*);
};


#endif
