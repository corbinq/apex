/*
	Corbin's notes on setOptions :
	
	GLOBAL OPTIONS
	
	setOptions is for setting global options, which are set at 
	runtime and used throughout (potentially) all source files
*/

#include "setOptions.hpp"

using namespace std;


namespace global_opts
{
	
	// GENERAL OPTIONS
		string out_prefix = "";
		bool low_mem = false;
	
	// INPUT OPTIONS 
		
		// GENOTYPE OPTIONS
		// bool set_from_
		int exclude_missing = 0;
		int set_missing_to_reference = 1;
		double minimum_maf = 0;
		int minimum_mac = 1;
		bool use_dosages = false;
		
		// COVARIATE OPTIONS 
		bool filter_covariates = false;
		// use_covariates;
		
		// SAMPLE SUBSETTING
		bool filter_iids = false;
		// include_iids;
		// exclude_iids;
		
		bool trim_gene_ids = false;
		
	// ANALYSIS OPTIONS 

		double RSQ_PRUNE;
		double RSQ_BUDDY = 1.00;

		// TESTING OPTIONS
		bool het_use_hom = true;
		bool het_use_het = true;
		bool het_use_acat = true;
		
		bool step_marginal = false;

		// VARIANT MERGE OPTIONS
		bool biallelic_only = true;

		// LMM OPTIONS
		bool use_grm = false;
		bool ml_not_reml = false;
		
		// ANALYSIS MODE
		bool IVW_H1_SIGMA = false;
		bool conditional_analysis = false;
		bool trans_eqtl_mode = false;
		bool backward_step = false;
		
		// CIS-eQTL OPTIONS
		double LM_ALPHA = 0.05;
		int cis_window_bp = 1000000;
		bool cis_window_gene_body = false;
			
		// LD OPTIONS
		int ld_window_bp = 1000000;
		
	// VARIANT MATCHING OPTIONS
	
		double freq_tol = 0.05;
		bool try_match_ambiguous_snv = true;
		
	// GENE SUBSETTING OPTIONS
	
		bool filter_genes;
		vector<string> target_genes;
		
	// GLOBAL VALUES
	
		unordered_map<string,int> i_chrom_map({
			{"1" , 1},{"2" , 2},{"3" , 3},{"4" , 4},
			{"5" , 5},{"6" , 6},{"7" , 7},{"8" , 8},
			{"9" , 9},{"10" , 10},{"11" , 11},{"12" , 12},
			{"13" , 13},{"14" , 14},{"15" , 15},{"16" , 16},
			{"17" , 17},{"18" , 18},{"19" , 19},{"20" , 20},
			{"21" , 21},{"22" , 22},{"X" , 23},{"Y" , 24},
			{"M" , 25},{"MT" , 25}, 
			{"chr1" , 1},{"chr2" , 2},{"chr3" , 3},{"chr4" , 4},
			{"chr5" , 5},{"chr6" , 6},{"chr7" , 7},{"chr8" , 8},
			{"chr9" , 9},{"chr10" , 10},{"chr11" , 11},{"chr12" , 12},
			{"chr13" , 13},{"chr14" , 14},{"chr15" , 15},{"chr16" , 16},
			{"chr17" , 17},{"chr18" , 18},{"chr19" , 19},{"chr20" , 20},
			{"chr21" , 21},{"chr22" , 22},{"chrX" , 23},{"chrY" , 24},
			{"chrM" , 25},{"chrMT" , 25}
		});
}

bool global_opts::process_global_opts( const string& pfx, const bool& low_memory, const double& rsq_buddy, const double& rsq, const double& pthresh, const int& window, const vector<string>& tg, const bool& ivw_mode, const bool& use_ds, const bool& trim, const bool& backward, const bool& h_hom, const bool& h_het, const bool& h_acat, const bool& step_marg ){
	out_prefix = pfx;
	low_mem = low_memory;
	cis_window_bp = window;
	ld_window_bp = 2*window;
	RSQ_BUDDY = rsq_buddy;
	RSQ_PRUNE = rsq;
	LM_ALPHA = pthresh;
	filter_genes = (tg.size() > 0);
	target_genes = tg;
	IVW_H1_SIGMA = ivw_mode;
	use_dosages = use_ds;
	trim_gene_ids = trim;
	backward_step = backward;
	het_use_hom = h_hom;
	het_use_het = h_het;
	het_use_acat = h_acat;
	step_marginal = step_marg;
	
	return true;
};



