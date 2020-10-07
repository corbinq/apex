# 10 "yax/src/setOptions.hpp.c"
#ifndef SETOPTIONS_HPP
#define SETOPTIONS_HPP 

#include <iostream>
#include <vector>
#include <unordered_map>

namespace global_opts
{

  extern std::string out_prefix;
  extern bool low_mem;




  extern int exclude_missing;
  extern int set_missing_to_reference;
  extern double minimum_maf;
  extern int minimum_mac;
  extern bool use_dosages;


  extern bool filter_covariates;
  extern std::vector<std::string> use_covariates;


  extern bool filter_iids;
  extern std::vector<std::string> include_iids;
  extern std::vector<std::string> exclude_iids;


  extern bool trim_gene_ids;



  extern double RSQ_PRUNE;
  extern double RSQ_BUDDY;


  extern bool het_use_hom;
  extern bool het_use_het;
  extern bool het_use_acat;

  extern bool step_marginal;


  extern bool biallelic_only;


  extern bool use_grm;
  extern bool ml_not_reml;


  extern bool IVW_H1_SIGMA;
  extern bool conditional_analysis;
  extern bool trans_eqtl_mode;
  extern bool backward_step;


  extern double LM_ALPHA;
  extern int cis_window_bp;
  extern bool cis_window_gene_body;


  extern int ld_window_bp;



  extern bool filter_genes;
  extern std::vector<std::string> target_genes;



  extern double freq_tol;
  extern bool try_match_ambiguous_snv;



  extern std::unordered_map<std::string,int> i_chrom_map;



  bool process_global_opts(const std::string& pfx, const bool& low_memory, const double& rsq_buddy, const double& rsq, const double& pthresh, const int& window, const std::vector<std::string>& tg, const bool& ivw_mode, const bool& use_ds, const bool& trim, const bool& backward, const bool& h_hom, const bool& h_het, const bool& h_acat, const bool& step_marg);
}

# endif
