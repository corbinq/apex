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


#ifndef CISMAPPING_HPP
#define CISMAPPING_HPP

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


void run_cis_QTL_analysis(bcf_srs_t*& sr, bcf_hdr_t*& hdr, genotype_data& g_data, table& c_data, bed_data& e_data, block_intervals& bm, const bool& rknorm_y, const bool& rknorm_r, const bool& make_sumstat, const bool& make_long, const bool& just_long);

void run_cis_QTL_analysis_LMM(bcf_srs_t*& sr, bcf_hdr_t*& hdr,genotype_data& g_data, table& c_data, bed_data& e_data, Eigen::SparseMatrix<double>& GRM, const std::vector<int>& relateds, block_intervals& bm, const bool& rknorm_y, const bool& rknorm_r, const bool& make_sumstat, const bool& make_long, const bool& just_long);

void run_cis_QTL_analysis_eLMM(const int& n_fac, const int& n_fac_fe, bcf_srs_t*& sr, bcf_hdr_t*& hdr,genotype_data& g_data, table& c_data, bed_data& e_data, block_intervals& bm, const bool& rknorm_y, const bool& rknorm_r, const bool& make_sumstat, const bool& make_long, const bool& just_long, Eigen::MatrixXd& Y_epc);

void run_cis_QTL_analysis_eFE (const int& n_fac, bcf_srs_t*& sr, bcf_hdr_t*& hdr,genotype_data& g_data, table& c_data, bed_data& e_data, block_intervals& bm, const bool& rknorm_y, const bool& rknorm_r, const bool& make_sumstat, const bool& make_long, const bool& just_long, Eigen::MatrixXd& Y_epc );

void scan_signals(bcf_srs_t*& sr, bcf_hdr_t*& hdr,genotype_data& g_data, table& c_data, bed_data& e_data, block_intervals& bm, const bool& rknorm_y, const bool& rknorm_r);

void save_fa_covariates(const int& n_fac, genotype_data& g_data, table& c_data, bed_data& e_data, const bool& rknorm_y, const bool& rknorm_r);

#endif

