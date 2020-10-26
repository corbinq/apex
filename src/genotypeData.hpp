/*  
    Copyright (C) 2020 
    Author: Corbin Quick <qcorbin@hsph.harvard.edu>

    This file is part of YAX.

    YAX is distributed "AS IS" in the hope that it will be 
    useful, but WITHOUT ANY WARRANTY; without even the implied 
    warranty of MERCHANTABILITY, NONINFRINGEMENT, or FITNESS 
    FOR A PARTICULAR PURPOSE.

    The above copyright notice and this permission notice shall 
    be included in all copies or substantial portions of YAX.
*/


#ifndef GENOTYPEDATA_HPP
#define GENOTYPEDATA_HPP

#include <vector>
#include <unordered_map>

#include "setOptions.hpp"
#include "mapID.hpp"
#include "htsWrappers.hpp"
#include "dataParser.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>


void read_sparse_GRM(const std::string& filename, Eigen::SparseMatrix<double>& GRM, const std::vector<std::string>& kp_ids, const double& r_scale, const int& r_col, std::vector<int>& related);

void read_dense_GRM(const std::string&, Eigen::MatrixXd&, std::vector<std::string>&);

/*
class rel_blocks
{
	public:
		void addPair(const int& i, const int&j){
			
		};
		std::
	
}
*/

class sparse_gt
{
	private:
		std::vector<std::pair<int, int>> entries;
		double sum_gt;
		int n_samples;
		int n_missing;
		int last;
	public:
		sparse_gt(const int& ns) : sum_gt(0.00), n_samples(ns), n_missing(0), last(-1) {};
		int size(){
			return entries.size();
		};
		void set_gt(const int& n, const int& gt){
			if( n > last ){
				if( gt < 0 ){
					entries.push_back(std::make_pair(n,-1));
					n_missing++;
				}else if( gt > 0 ){
					entries.push_back(std::make_pair(n,gt));
					sum_gt += gt;
				}
				last = n;
			}else{
				std::cerr << "Fatal: Tried to add gt out of order!\n";
				exit(2);
			}
		};
		void flip(const int& n){
			int j = 0, i = 0;
			std::vector<std::pair<int, int>> flipped;
			while(i < n){
				if( j < entries.size() ){
					while( i < entries[j].first && i < n ){
						flipped.push_back( std::make_pair(i, 2) );
						i++;
					}
					if( i == entries[j].first ){
						if( entries[j].second < 0 ){
							flipped.push_back(entries[j]);
						}else if( entries[j].second == 1 ){
							flipped.push_back(entries[j]);
						}else if(entries[j].second != 2){
							std::cerr <<"\n\n" << entries[j].second << "\n\n";
							std::cerr << "entries[j].second != 2\n";
							exit(1);
						}
						i++;
						j++;
					}else if( j > i ){
						std::cerr << "This should never happen!\n";
						exit(2);
					}
				}else{
					flipped.push_back( std::make_pair(i, 2) );
					i++;
				}
			}
			sum_gt = 2.00 * ( (double) (n_samples - n_missing) ) - sum_gt;
			entries = flipped;
		};
		void add_gt_sparsemat(Eigen::SparseMatrix<double>& smat, const int& col_n){
			smat.startVec(col_n);
			double missing_val = sum_gt/( (double) (n_samples - n_missing) );
			for( const auto& x : entries){
				if( x.second > 0 ){
					smat.insertBack(x.first, col_n) = x.second;
				}else{
					smat.insertBack(x.first, col_n) = missing_val;
				}
			}
		};
};


class sparse_ds
{
	public:
		std::vector<float> dosages;
		int last;
		double sum_dos;
		int n_missing;
		int n_samples;
		
		sparse_ds(const int& N) : sum_dos(0.00), n_missing(0), n_samples(N), last(-1) {dosages.resize(N);};
		
		void set_ds(const int& n, const float& ds){
			if( n > last ){
				if( ds < 0 ){
					n_missing++;
				}else{
					sum_dos += ds;
				}
				dosages[n] = ds;
				last = n;
			}else{
				std::cerr << "Fatal: Tried to add ds out of order!\n";
				exit(2);
			}
		};
		void flip(){
			for( float& ds: dosages ){
				if(ds >= 0){
					ds = 2.00 - ds;
				}
			}
			sum_dos = 2.00 * ( (double) (n_samples - n_missing) ) - sum_dos;
		};
		void add_ds_sparsemat(Eigen::SparseMatrix<double>& smat, const int& col_n){
			double missing_val = sum_dos/( (double) (n_samples - n_missing) );
			smat.startVec(col_n);
			for(int i = 0; i < dosages.size(); i++){
				if( dosages[i] < 0.00 ){
					if( missing_val > global_opts::dosage_thresh ){
						smat.insertBack(i, col_n) = missing_val;
					}
				}else if( dosages[i] > global_opts::dosage_thresh ){
					smat.insertBack(i, col_n) = (double) dosages[i];
				}
			}
		};
};


class genotype_data
{
	public:
		genotype_data(): chunk_size(1000), nz_frac(0.1), n_variants(0), n_samples(0) {};
	
		int chunk_size;
		double nz_frac;

		int n_samples;
		int n_variants;
		
		id_map ids;

		std::vector<std::string> chr;
		std::vector<int> pos;
		//std::vector<std::string> rsid; 
		std::vector<std::string> ref;
		std::vector<std::string> alt;
		
		std::vector<std::pair<int,int>> ld_index_range;
		
		std::vector<double> mean;
		std::vector<double> var;
		
		std::vector<bool> flipped;
		std::vector<bool> keep;
		
		int geno_start;
		int geno_size;
		
		Eigen::SparseMatrix<double> genotypes;
		Eigen::MatrixXd dense_genotypes;
		
		void read_bcf_variants(bcf_srs_t*&, bcf_hdr_t*&, int&, bool store_geno = true, bool scan_geno = true);
		
		int get_ld_index();
		inline bool process_bcf_variant(bcf1_t*&, bcf_hdr_t*&, bool store_geno = true, bool scan_geno = true);
		void read_bcf_header(bcf_hdr_t*);
		void freeze_genotypes();
		void resize_genotypes(const int& nv){
			genotypes.resize(ids.keep.size(), nv);
			genotypes.reserve( (int) (nz_frac * genotypes.rows() * genotypes.cols() ) );
		};
		void initialize_genotypes(const int& nv){
			resize_genotypes(nv);
		};
		void print_genotypes(){
			std::cout << genotypes << "\n";
		};
		
		void clear();
		void clear_genotypes();
		void read_genotypes(bcf_srs_t*&, bcf_hdr_t*&, const int&, const int&);
		
		bool record_matches(bcf1_t*& rec, bcf_hdr_t*& hdr, const int& i);
		
	private:
		inline bool add_bcf_genotypes(int*&, const int&, double&, double&, bool&, const bool);
		inline bool add_bcf_dosages(float*& ds_rec, const int& col_n, double& mean_, double& var_, bool& flip_, const bool store_geno);
		
		void checkResize(const int& col_n){
			if( col_n >= genotypes.cols() ){
				genotypes.conservativeResize(ids.keep.size(), col_n + chunk_size);
			}
		};
		
};


#endif

