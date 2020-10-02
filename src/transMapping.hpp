#ifndef TRANSMAPPING_HPP
#define TRANSMAPPING_HPP

#include "setOptions.hpp"
#include "readBED.hpp"
#include "setRegions.hpp"
#include "readTable.hpp"
#include "genotypeData.hpp"
#include "htsWrappers.hpp"
#include "fitUtils.hpp"
#include "miscUtils.hpp"
#include "mathStats.hpp"

#include "cisMapping.hpp"



class theta_data
{
	private:
		std::vector<std::string> chr;
		std::vector<int> start;
		std::vector<int> end;
		std::vector<std::string> gene_id;
		std::vector<double> sigma2;
		std::vector<double> tau2;
		std::vector<double> phi;
		
		std::unordered_map<std::string, int> gene_map; 
		
	public:
		theta_data(const std::string& prefix, const std::string& region = ""){
			std::vector<double> tau2;
			
			data_parser dp;
			
			dp.add_field(chr, 0);
			dp.add_field(start, 1);
			dp.add_field(end, 2);
			dp.add_field(gene_id, 3);
			dp.add_field(sigma2, 4);
			dp.add_field(tau2, 5);
			dp.add_field(phi, 6);
			
			dp.parse_file(prefix + ".theta.gz", region);
			
			for( int i = 0; i < gene_id.size(); i++){
				gene_map[gene_id[i]] = i;
			}
			return;
		};
		void getTheta(const int& i, const std::string& gene, double& sigma2_, double& tau2_, double& phi_){
			int j;
			if( gene_id[i] == gene ){
				j = i;
			}else{
				j = gene_map[gene];
			}
			sigma2_ = sigma2[j];
			tau2_ = tau2[j];
			phi_ = phi[j];
			return;
		};
};

void run_trans_eQTL_analysis(bcf_srs_t*& sr, bcf_hdr_t*& hdr, genotype_data& g_data, table& c_data, bed_data& e_data, const bool& rknorm_y, const bool& rknorm_r, const bool& make_sumstat, const bool& make_long, const int& chunk_size);

void run_trans_eQTL_analysis_LMM(bcf_srs_t*& sr, bcf_hdr_t*& hdr, genotype_data& g_data, table& c_data, bed_data& e_data, Eigen::SparseMatrix<double>& GRM, const std::vector<int>& relateds, const bool& rknorm_y, const bool& rknorm_r, const bool& make_sumstat, const bool& make_long, const int& chunk_size);

#endif

