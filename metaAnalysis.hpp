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

using namespace std;

class cis_sumstat_data
{
	public:
		string file_prefix;
		string region;
	
		lindex ln;
	
		vector<string> chr; 
		vector<int> start; 
		vector<int> end;
		vector<string> gene_id;
		
		vector<double> egene_pval;
		
		vector<int> NS;
		vector<int> NC;
		vector<double> SD;
		
		vector<int> S_CIS;
		vector<int> N_CIS; 
		
		vector<Eigen::VectorXd> score;
		
		int n_genes;

		cis_sumstat_data(const string&,const string& reg = "");
		void open(const string&,const string& reg = "");
};

class vcov_meta_data
{
	public:
		vector<vcov_data> vc;
		
		vector<string> chr;
		vector<int> pos;
		vector<string> ref;
		vector<string> alt;
		
		vector<double> mac;
		vector<double> var;
		
		vector<vector<double>> var_perStudy;
		
		// indexes for shared variants within each study
			// si[i][j] = the index of the jth shared snp
			//                     within the list of snps for study i 
		vector<vector<int>> si; 
		
		// one vector per study, indicating whether flipped relative to study #1
		vector<vector<bool>> sflipped;
		
		vcov_meta_data() {};
		vcov_meta_data(const vector<string>& v_fn, const string& reg = ""){
			open(v_fn, reg);
		};
		
		void open(const vector<string>& v_fn, const string& reg = "");
		
		Eigen::MatrixXd getGtG(const vector<int>&, const vector<double>& w = vector<double>(0), const vector<int>& sl = vector<int>(0) );
		Eigen::MatrixXd getV(const vector<int>&, const vector<double>& w = vector<double>(0), const vector<int>& sl = vector<int>(0) );
		
		Eigen::MatrixXd getGtG(const vector<int>&,const vector<int>&, const vector<double>& w = vector<double>(0), const vector<int>& sl = vector<int>(0) );
		Eigen::MatrixXd getV(const vector<int>&,const vector<int>&, const vector<double>& w = vector<double>(0), const vector<int>& sl = vector<int>(0) );
		
		Eigen::MatrixXd getGtG_perStudy(const int&, const vector<int>&);
		Eigen::MatrixXd getV_perStudy(const int&, const vector<int>&);
		
		Eigen::MatrixXd getGtG_perStudy(const int&, const vector<int>&,const vector<int>&);
		Eigen::MatrixXd getV_perStudy(const int&, const vector<int>&,const vector<int>&);
		
	private:
		void build_index();
};


class cis_meta_data
{
	public:
		vector<cis_sumstat_data> ss;
		
		vcov_meta_data vc;
		
		vector<string> chr; 
		vector<int> start; 
		vector<int> end;
		vector<string> gene_id;
		
		vector<vector<int>> study_list;
		
		vector<double> egene_pval;

		vector<double> N;
		vector<double> DF;
		vector<double> SD;
		vector<double> ADJ;
		
		vector<int> S_CIS;
		vector<int> N_CIS;
		
		vector<vector<double>> ivw;
		vector<vector<double>> SD_perStudy;
		vector<vector<double>> SSR_perStudy;
		vector<vector<double>> DF_perStudy;
		
		vector<Eigen::VectorXd> score;
		vector<Eigen::VectorXd> var_score;
		
		vector<vector<Eigen::VectorXd>> score_perStudy;
		vector<vector<Eigen::VectorXd>> var_score_perStudy;
		
		cis_meta_data() {};
		cis_meta_data(const vector<string>& v_fn, const string& reg = "" ){
			vc.open(v_fn, reg);
			open(v_fn, reg);
			merge(vc.si, vc.sflipped);
		};
		
		void merge_intersection(const vector<vector<int>>&, const vector<vector<bool>>&);
		void merge(const vector<vector<int>>&, const vector<vector<bool>>&);
		void meta_analyze();
		void conditional_analysis(const int&, ostream&, ostream&);
		void conditional_analysis();
		void conditional_analysis_het(const int&, ostream&);
		void conditional_analysis_het();
		
		void add_study(const string& fn, const string& reg = ""){
			ss.push_back( cis_sumstat_data(fn, reg) );
		};
		void open(const vector<string>& v_fn, const string& reg = ""){
			for( const string& fn : v_fn ){ add_study(fn,reg); }
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
				// cout << SD_perStudy[g][s] << ", " << vc.var_perStudy[S_CIS[g] + i][s] << "\n";
				dv_val += diagV_perStudy(g,i,s) * ivw[g][s];
			}
			// cout << dv_val << "\n";
			return dv_val;
		}
		
		Eigen::VectorXd dV_perStudy(const int& g, const int& s){
			Eigen::VectorXd out(N_CIS[g]);
			for( int i = 0; i < N_CIS[g]; i++){
				out(i) = diagV_perStudy(g, i, s);
			}
			return out;
		}
		
		vector<double> ivw_weight_H1(const int& g, const int& i){
			vector<double> ivw_H1(ss.size());
			for( const int& s : study_list[g] ){
				double sc = score_perStudy[g][s](i);
				ivw_H1[s] = (DF_perStudy[g][s]-1)/(SSR_perStudy[g][s] - sc*sc/diagV_perStudy(g,i,s) );
			}
			return ivw_H1;
		}
		
		const int& pos(const int& i){ return vc.pos[i];};
		const string& ref(const int& i){ return vc.ref[i];};
		const string& alt(const int& i){ return vc.alt[i];};
};


class vcov_getter
{
	private:
		vector<double>& ivw;
		vcov_meta_data& vc;
		int s_var;
		int n_var;
		
		vector<int> use_studies;
		
		vector<int> R_CIS;
		vector<int> add(vector<int> inp, const int& s){
			for( int &i : inp ){i += s;}
			return inp;
		}
	public:
		vcov_getter(vcov_meta_data& vc_, vector<double>& ivw_, int s_var_, int n_var_) : 
			vc(vc_), ivw(ivw_), s_var(s_var_), n_var(n_var_)
		{
			for(int i = 0; i < n_var; i++){
				R_CIS.push_back(s_var + i);
			}
			for(int s = 0; s < vc.si.size(); s++ ){
				use_studies.push_back(s);
			}
		};
		vcov_getter(vcov_meta_data& vc_, vector<double>& ivw_, int s_var_, int n_var_, vector<int> use_studies_) : 
			vc(vc_), ivw(ivw_), s_var(s_var_), n_var(n_var_), use_studies(use_studies_)
		{
			for(int i = 0; i < n_var; i++){
				R_CIS.push_back(s_var + i);
			}
		};
		Eigen::VectorXd Covar(const int& i){ return vc.getV(R_CIS, vector<int>(1,s_var + i), ivw, use_studies);};
		Eigen::MatrixXd Var(const vector<int>& wh){ return vc.getV(add(wh, s_var), ivw, use_studies);};
		
		Eigen::VectorXd Covar_perStudy(const int& s, const int& i){ return vc.getV_perStudy(use_studies[s], R_CIS, vector<int>(1,s_var + i));};
		Eigen::MatrixXd Var_perStudy(const int& s, const vector<int>& wh){ return vc.getV_perStudy(use_studies[s], add(wh, s_var) );};
};


#endif

