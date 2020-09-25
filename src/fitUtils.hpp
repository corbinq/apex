/*
	Centering, scaling, residualizing, and transformations
*/

#ifndef FITUTILS_HPP
#define FITUTILS_HPP

#include <math.h>
#include <vector>
#include <functional>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "mathStats.hpp"


void appendInterceptColumn(Eigen::MatrixXd&);
Eigen::MatrixXd get_half_hat_matrix(const Eigen::MatrixXd&);
Eigen::MatrixXd resid_from_half_hat( Eigen::MatrixXd, const Eigen::MatrixXd&);

void scale_and_center(Eigen::MatrixXd&, std::vector<double>&);
void scale_and_center(Eigen::MatrixXd&);
void rank_normalize (Eigen::MatrixXd&);

#endif

