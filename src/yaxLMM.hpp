/*  
    Copyright (C) 2020 
    Author: Corbin Quick <qcorbin@hsph.harvard.edu>

    This file is a part of YAX.

    YAX is distributed "AS IS" in the hope that it will be 
    useful, but WITHOUT ANY WARRANTY; without even the implied 
    warranty of MERCHANTABILITY, NON-INFRINGEMENT, or FITNESS 
    FOR A PARTICULAR PURPOSE.

    The above copyright notice and disclaimer of warranty must 
    be included in all copies or substantial portions of YAX.
*/


#ifndef YAXLMM_HPP
#define YAXLMM_HPP

#include "setOptions.hpp"
#include "readBED.hpp"
#include "setRegions.hpp"
#include "readTable.hpp"
#include "genotypeData.hpp"
#include "htsWrappers.hpp"
#include "fitUtils.hpp"
#include "fitModels.hpp"
#include "miscUtils.hpp"
#include "mathStats.hpp"

void fit_LMM_null_models(table& c_data, bed_data& e_data, Eigen::SparseMatrix<double>& GRM, const std::vector<int>& relateds, const bool& rknorm_y, const bool& rknorm_r);

void calculate_V_anchor_points( Eigen::MatrixXd& V_mat, genotype_data& g_data, const Eigen::MatrixXd& C, const std::vector<double>& hsq_vals, const Eigen::MatrixXd& CtC, const Eigen::MatrixXd& CtC_i, const Eigen::MatrixXd& QtG, Eigen::MatrixXd& QtC, const Eigen::VectorXd Q_lambda );

void calculate_V_anchor_points_low_rank( Eigen::MatrixXd& V_mat, genotype_data& g_data, const Eigen::MatrixXd& C, const std::vector<double>& hsq_vals, const Eigen::MatrixXd& CtC, const Eigen::MatrixXd& CtC_i, const Eigen::MatrixXd& QtG, Eigen::MatrixXd& QtC, const Eigen::VectorXd Q_lambda );

#endif
