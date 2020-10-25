/*  
    Copyright (C) 2020 
    Author: Corbin Quick <qcorbin@hsph.harvard.edu>

    This file is part of YAX.

    YAX is distributed "AS IS" in the hope that it will be 
    useful, but WITHOUT ANY WARRANTY; without even the implied 
    warranty of MERCHANTABILITY, NONINFRINGEMENT, or FITNESS 
    FOR A PARTICULAR PURPOSE.

    The above copyright notice and this permission notice shall 
    be included in all copies or substantial portions of YAX.
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
Eigen::MatrixXd resid_from_half_hat(const Eigen::MatrixXd&, const Eigen::MatrixXd&);

Eigen::VectorXd resid_vec_from_half_hat(const Eigen::VectorXd& Y, const Eigen::MatrixXd&);
void make_half_hat_matrix(Eigen::MatrixXd& X);
void make_resid_from_half_hat( Eigen::MatrixXd& Y, const Eigen::MatrixXd& C );


void scale_and_center(Eigen::MatrixXd&, std::vector<double>&);
void scale_and_center(Eigen::MatrixXd&);
void rank_normalize (Eigen::MatrixXd&);

#endif

