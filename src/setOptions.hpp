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


/*
	setOptions is for setting global options, which are set at 
	runtime and used throughout multiple source files
*/

#ifndef SETOPTIONS_HPP
#define SETOPTIONS_HPP

#include <iostream>
#include <vector>
#include <unordered_map>

namespace global_opts
{
	// GENERAL OPTIONS
		extern std::string out_prefix;
		extern bool low_mem;
    extern bool write_logp;
    extern bool legacy_vcov;

	// INPUT OPTIONS

		extern std::string global_region;

		// GENOTYPE OPTIONS
		extern int exclude_missing;
		extern double dosage_thresh;
		extern double minimum_maf;
		extern int minimum_mac;
		extern bool use_dosages;

		// COVARIATE OPTIONS
		extern bool filter_covariates;
		extern std::vector<std::string> use_covariates;

		// SAMPLE SUBSETTING
		extern bool filter_iids;
		extern std::vector<std::string> include_iids;
		extern std::vector<std::string> exclude_iids;

		// GENE OPTIONS
		extern bool trim_gene_ids;

	// ANALYSIS OPTIONS

		extern bool write_resid_mat;

		extern int n_fa_iter;
		extern double fa_p;
		extern double fa_tau;

		// stepwise options
		extern int max_signals;
		extern int max_steps;

		extern double exp_weight_val;

		// "Sloppy" covariate adjustment
		extern bool sloppy_covar;

		extern double RSQ_PRUNE;
		extern double RSQ_BUDDY;

		// TESTING OPTIONS
		extern bool het_use_hom;
		extern bool het_use_het;
		extern bool het_use_acat;

		extern bool step_marginal;

		// VARIANT MERGE OPTIONS
		extern bool biallelic_only;

		// LMM OPTIONS
		extern bool use_grm;
		extern bool ml_not_reml;
		extern bool write_v_anchors;

		// ANALYSIS MODE
		extern char meta_weight_method;
		extern bool conditional_analysis;
		extern bool trans_eqtl_mode;
		extern double backward_thresh;

		// CIS-QTL OPTIONS
		extern double LM_ALPHA;
		extern int cis_window_bp;
		extern bool cis_window_gene_body;

		// LD OPTIONS
		extern int ld_window_bp;

	// GENE SUBSETTING OPTIONS

		extern bool filter_genes;
		extern std::vector<std::string> target_genes;

	// VARIANT MATCHING OPTIONS

		extern double freq_tol;
		extern bool try_match_ambiguous_snv;

	// GLOBAL VARIABLES

		extern std::unordered_map<std::string,int> i_chrom_map;

	// PROCESS OPTIONS

		bool process_global_opts(const std::string& pfx, const bool& low_memory, const double& rsq_buddy, const double& rsq, const double& pthresh, const int& window, const std::vector<std::string>& tg, const char& ivw_mode, const bool& use_ds, const bool& trim, const double& backward, const bool& h_hom, const bool& h_het, const bool& h_acat, const bool& step_marg, const bool& print_logp);

		bool set_lmm_options(const bool& wap);

		void set_max_signals(const int& ms);

		void use_sloppy_covar();

		void set_exp_weight(const double&);

		void save_residuals(const bool&);

		void set_factor_par(const int&, const double&, const double&);

		void set_global_region(const std::string& );
    void set_legacy_vcov(const bool& b);
    void reset();
}
# endif

