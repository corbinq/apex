#ifndef METAANALYSIS_HPP
#define METAANALYSIS_HPP

#include <vector>
#include <fstream> 

#include "setOptions.hpp"
#include "processVCOV.hpp"
#include "miscUtils.hpp"
#include "fitModels.hpp"
#include "dataParser.hpp"

#include <Eigen/Sparse>
#include <Eigen/Dense>


class cis_sumstat_data
{
	public:
		std::string file_prefix;
		std::string region;
	
		lindex ln;
	
		std::vector<std::string> chr; 
		std::vector<int> start; 
		std::vector<int> end;
		std::vector<std::string> gene_id;
		
		std::vector<double> egene_pval;
		
		std::vector<int> NS;
		std::vector<int> NC;
		std::vector<double> SD;
		
		std::vector<int> S_CIS;
		std::vector<int> N_CIS; 
		
		std::vector<Eigen::VectorXd> score;
		
		int n_genes;

		cis_sumstat_data(const std::string&,const std::string& reg = "");
		void open(const std::string&,const std::string& reg = "");
};

class vcov_meta_data
{
	public:
		std::vector<vcov_data> vc;
		
		std::vector<std::string> chr;
		std::vector<int> pos;
		std::vector<std::string> ref;
		std::vector<std::string> alt;
		
		std::vector<double> mac;
		std::vector<double> var;
		
		std::vector<std::vector<double>> var_perStudy;
		
		// indexes for shared variants within each study
			// si[i][j] = the index of the jth shared snp
			//                     within the list of snps for study i 
		std::vector<std::vector<int>> si; 
		
		// one vector per study, indicating whether flipped relative to study #1
		std::vector<std::vector<bool>> sflipped;
		
		vcov_meta_data() {};
		vcov_meta_data(const std::vector<std::string>& v_fn, const std::string& reg = ""){
			open(v_fn, reg);
		};
		
		void open(const std::vector<std::string>& v_fn, const std::string& reg = "");
		
		Eigen::MatrixXd getGtG(const std::vector<int>&, const std::vector<double>& w = std::vector<double>(0), const std::vector<int>& sl = std::vector<int>(0) );
		Eigen::MatrixXd getV(const std::vector<int>&, const std::vector<double>& w = std::vector<double>(0), const std::vector<int>& sl = std::vector<int>(0) );
		
		Eigen::MatrixXd getGtG(const std::vector<int>&,const std::vector<int>&, const std::vector<double>& w = std::vector<double>(0), const std::vector<int>& sl = std::vector<int>(0) );
		Eigen::MatrixXd getV(const std::vector<int>&,const std::vector<int>&, const std::vector<double>& w = std::vector<double>(0), const std::vector<int>& sl = std::vector<int>(0) );
		
		Eigen::MatrixXd getGtG_perStudy(const int&, const std::vector<int>&);
		Eigen::MatrixXd getV_perStudy(const int&, const std::vector<int>&);
		
		Eigen::MatrixXd getGtG_perStudy(const int&, const std::vector<int>&,const std::vector<int>&);
		Eigen::MatrixXd getV_perStudy(const int&, const std::vector<int>&,const std::vector<int>&);
		
	private:
		void build_index();
};


class cis_meta_data
{
	public:
		std::vector<cis_sumstat_data> ss;
		
		vcov_meta_data vc;
		
		std::vector<std::string> chr; 
		std::vector<int> start; 
		std::vector<int> end;
		std::vector<std::string> gene_id;
		
		std::vector<std::vector<int>> study_list;
		
		std::vector<double> egene_pval;

		std::vector<double> N;
		std::vector<double> DF;
		std::vector<double> SD;
		std::vector<double> ADJ;
		
		std::vector<int> S_CIS;
		std::vector<int> N_CIS;
		
		std::vector<std::vector<double>> ivw;
		std::vector<std::vector<double>> SD_perStudy;
		std::vector<std::vector<double>> SSR_perStudy;
		std::vector<std::vector<double>> DF_perStudy;
		
		std::vector<Eigen::VectorXd> score;
		std::vector<Eigen::VectorXd> var_score;
		
		std::vector<std::vector<Eigen::VectorXd>> score_perStudy;
		std::vector<std::vector<Eigen::VectorXd>> var_score_perStudy;
		
		cis_meta_data() {};
		cis_meta_data(const std::vector<std::string>& v_fn, const std::string& reg = "" ){
			vc.open(v_fn, reg);
			open(v_fn, reg);
			merge(vc.si, vc.sflipped);
		};
		
		void merge_intersection(const std::vector<std::vector<int>>&, const std::vector<std::vector<bool>>&);
		void merge(const std::vector<std::vector<int>>&, const std::vector<std::vector<bool>>&);
		void meta_analyze();
		
		void conditional_analysis(const int&, std::ostream&, std::ostream&);
		void conditional_analysis();
		void conditional_analysis_het(const int&, std::ostream&);
		void conditional_analysis_het();
		
		void get_vcov_gene(const int& gene_index, const bool& centered = true);
		
		void add_study(const std::string& fn, const std::string& reg = ""){
			ss.push_back( cis_sumstat_data(fn, reg) );
		};
		void open(const std::vector<std::string>& v_fn, const std::string& reg = ""){
			for( const std::string& fn : v_fn ){ add_study(fn,reg); }
		};
		
	private:
		// const double& diagV(const int& i){return vc.var[i];};
		// const int& pos(const int& i){ return vc.pos[i];};
		
		double& diagV_perStudy(const int& g, const int& i, const int& s){
			return vc.var_perStudy[S_CIS[g] + i][s];
		}
		
		double diagV(const int& g, const int& i){
			double dv_val = 0;
			for( const int& s : study_list[g] ){
				// std::cout << SD_perStudy[g][s] << ", " << vc.var_perStudy[S_CIS[g] + i][s] << "\n";
				dv_val += diagV_perStudy(g,i,s) * ivw[g][s];
			}
			// std::cout << dv_val << "\n";
			return dv_val;
		}
		
		Eigen::VectorXd dV_perStudy(const int& g, const int& s){
			Eigen::VectorXd out(N_CIS[g]);
			for( int i = 0; i < N_CIS[g]; i++){
				out(i) = diagV_perStudy(g, i, s);
			}
			return out;
		}
		
		std::vector<double> ivw_weight_H1(const int& g, const int& i){
			std::vector<double> ivw_H1(ss.size());
			for( const int& s : study_list[g] ){
				double sc = score_perStudy[g][s](i);
				ivw_H1[s] = (DF_perStudy[g][s]-1)/(SSR_perStudy[g][s] - sc*sc/diagV_perStudy(g,i,s) );
			}
			return ivw_H1;
		}
		
		const int& pos(const int& i){ return vc.pos[i];};
		const std::string& ref(const int& i){ return vc.ref[i];};
		const std::string& alt(const int& i){ return vc.alt[i];};
};


class vcov_getter
{
	private:
		std::vector<double>& ivw;
		vcov_meta_data& vc;
		int s_var;
		int n_var;
		
		std::vector<int> use_studies;
		
		std::vector<int> R_CIS;
		std::vector<int> add(std::vector<int> inp, const int& s){
			for( int &i : inp ){i += s;}
			return inp;
		}
	public:
		vcov_getter(vcov_meta_data& vc_, std::vector<double>& ivw_, int s_var_, int n_var_) : 
			vc(vc_), ivw(ivw_), s_var(s_var_), n_var(n_var_)
		{
			for(int i = 0; i < n_var; i++){
				R_CIS.push_back(s_var + i);
			}
			for(int s = 0; s < vc.si.size(); s++ ){
				use_studies.push_back(s);
			}
		};
		vcov_getter(vcov_meta_data& vc_, std::vector<double>& ivw_, int s_var_, int n_var_, std::vector<int> use_studies_) : 
			vc(vc_), ivw(ivw_), s_var(s_var_), n_var(n_var_), use_studies(use_studies_)
		{
			for(int i = 0; i < n_var; i++){
				R_CIS.push_back(s_var + i);
			}
		};
		Eigen::VectorXd Covar(const int& i){ return vc.getV(R_CIS, std::vector<int>(1,s_var + i), ivw, use_studies);};
		Eigen::MatrixXd Var(const std::vector<int>& wh){ return vc.getV(add(wh, s_var), ivw, use_studies);};
		Eigen::MatrixXd Var_uncentered(const std::vector<int>& wh){ return vc.getGtG(add(wh, s_var), ivw, use_studies);};
		
		Eigen::VectorXd Covar_perStudy(const int& s, const int& i){ return vc.getV_perStudy(use_studies[s], R_CIS, std::vector<int>(1,s_var + i));};
		Eigen::MatrixXd Var_perStudy(const int& s, const std::vector<int>& wh){ return vc.getV_perStudy(use_studies[s], add(wh, s_var) );};
};


#endif

