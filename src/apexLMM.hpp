/*  
    Copyright (C) 2020 
    Author: Corbin Quick <qcorbin@hsph.harvard.edu>

    This file is a part of APEX.

    APEX is distributed "AS IS" in the hope that it will be 
    useful, but WITHOUT ANY WARRANTY; without even the implied 
    warranty of MERCHANTABILITY, NON-INFRINGEMENT, or FITNESS 
    FOR A PARTICULAR PURPOSE.

    The above copyright notice and disclaimer of warranty must 
    be included in all copies or substantial portions of APEX.
*/


#ifndef APEXLMM_HPP
#define APEXLMM_HPP

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

void fit_LMM_null_models(table& c_data, bed_data& e_data, Eigen::SparseMatrix<double>& L, Eigen::VectorXd& GRM_lambda, const bool& rknorm_y, const bool& rknorm_r);

void fit_LMM_null_models_low_rank(const int& n_fac, table& c_data, bed_data& e_data, const bool& rknorm_y, const bool& rknorm_r, Eigen::MatrixXd& Y_epc );

void calculate_V_anchor_points( Eigen::MatrixXd& V_mat, genotype_data& g_data, const Eigen::MatrixXd& C, const std::vector<double>& hsq_vals, const Eigen::MatrixXd& CtC, const Eigen::MatrixXd& CtC_i, const Eigen::MatrixXd& QtG, Eigen::MatrixXd& QtC, const Eigen::VectorXd Q_lambda );

void calculate_V_anchor_points_low_rank( Eigen::MatrixXd& V_mat, genotype_data& g_data, const Eigen::MatrixXd& C, const std::vector<double>& hsq_vals, const Eigen::MatrixXd& CtC, const Eigen::MatrixXd& CtC_i, const Eigen::MatrixXd& QtG, Eigen::MatrixXd& QtC, const Eigen::VectorXd Q_lambda );

void save_V_anchor_points(bcf_srs_t*& sr, bcf_hdr_t*& hdr, genotype_data& g_data, table& c_data, Eigen::SparseMatrix<double>& L, Eigen::VectorXd& GRM_lambda);

void read_V_anchor_points( const std::string& anchor_path, Eigen::MatrixXd& V_mat, genotype_data& g_data, const std::vector<double>& hsq_vals);

#endif
