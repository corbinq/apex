/*
	Corbin's notes on setOptions :
	
	GLOBAL OPTIONS
	
	setOptions is for setting global options, which are set at 
	runtime and used throughout (potentially) all source files
*/

#ifndef SETOPTIONS_HPP
#define SETOPTIONS_HPP

#include <iostream>
#include <vector>
#include <unordered_map>

using namespace std;

namespace global_opts
{
	// GENERAL OPTIONS
		extern string out_prefix;
		extern bool low_mem;
	
	// INPUT OPTIONS
	
		// GENOTYPE OPTIONS 
		extern int exclude_missing;
		extern int set_missing_to_reference;
		extern double minimum_maf;
		extern int minimum_mac;
		extern bool use_dosages;
		
		// COVARIATE OPTIONS 
		extern bool filter_covariates;
		extern vector<string> use_covariates;
		
		// SAMPLE SUBSETTING
		extern bool filter_iids;
		extern vector<string> include_iids;
		extern vector<string> exclude_iids;
		
		// GENE OPTIONS
		extern bool trim_gene_ids;
		
	// ANALYSIS OPTIONS 

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
		
		// ANALYSIS MODE
		extern bool IVW_H1_SIGMA;
		extern bool conditional_analysis;
		extern bool trans_eqtl_mode;
		extern bool backward_step;
		
		// CIS-eQTL OPTIONS
		extern double LM_ALPHA;
		extern int cis_window_bp;
		extern bool cis_window_gene_body;
	
		// LD OPTIONS
		extern int ld_window_bp;
		
	// GENE SUBSETTING OPTIONS
	
		extern bool filter_genes;
		extern vector<string> target_genes;
		
	// VARIANT MATCHING OPTIONS
	
		extern double freq_tol;
		extern bool try_match_ambiguous_snv;
		
	// GLOBAL VARIABLES
		
		extern unordered_map<string,int> i_chrom_map;
		
	// PROCESS OPTIONS
    
		bool process_global_opts(const string& pfx, const bool& low_memory, const double& rsq_buddy, const double& rsq, const double& pthresh, const int& window, const vector<string>& tg, const bool& ivw_mode, const bool& use_ds, const bool& trim, const bool& backward, const bool& h_hom, const bool& h_het, const bool& h_acat, const bool& step_marg);
}

# endif

