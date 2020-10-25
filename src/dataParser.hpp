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


/*
This header provides a generic interface for parsing multiple,  
heterogeneous fields from htsFiles into multiple, heterogeneous 
data objects and formats (std::vector<T> and Eigen::Matrix). 

See readBed.cpp for an example. 	
*/

#ifndef DATAPARSER_HPP
#define DATAPARSER_HPP

#include <vector>

#include "htsWrappers.hpp"
#include "setOptions.hpp"


static int int0 = 0;

// class for filling Eigen matrix entries by reference from file
class matrix_filler
{
	public:
		matrix_filler(): n(-1), chunk_size(5000), unknown_dim(false) {};
		
		// Unknown matrix size (unknown number of input columns)
		void set(Eigen::MatrixXd* new_mat, const bool is_trans, int i_offset){
			mat = new_mat;
			start_at = i_offset;
			is_transposed = is_trans;
			n = 0;
			unknown_dim = true;
		}
		
		// Known matrix size; assign columns from interval range 
		void set(Eigen::MatrixXd* new_mat, const bool is_trans, const int& start, const int& n_fields){
			mat = new_mat;
			for( int i = start; i < start + n_fields; i++){
				idx.push_back(i);
			}
			is_transposed = is_trans;
			n = 0;
			if( n_fields != ( is_trans ? mat->rows() : mat->cols() ) ){
				std::cerr << "Fatal: matrix_filler error, mismatched dimensions\n";
				std::cerr << n_fields << " != (" << mat->cols() << ", " << mat->rows() << ").\n";
				abort();
			}
		}
		
		// Known matrix size; assign columns from index vector 
		void set(Eigen::MatrixXd* new_mat, const bool is_trans, const std::vector<int>& kp_index, int i_offset = 0){
			mat = new_mat;
			idx = kp_index;
			if( i_offset != 0 ){
				for(int j = 0; j < idx.size(); j++){
					idx[j] += i_offset;
				}
			}
			is_transposed = is_trans;
			n = 0;
			if( kp_index.size() != ( is_trans ? mat->rows() : mat->cols() ) ){
				std::cerr << "Fatal: matrix_filler error, mismatched dimensions\n";
				std::cerr << kp_index.size() << " != (" << mat->cols() << ", " << mat->rows() << ").\n";
				abort();
			}
		}
		
		void push_line(kstring_t& str, int*& offsets, int& n_fields){
			if( unknown_dim ){
				setSize(n_fields);
			}else{
				checkResize();
			}
			if( n >= 0 ){
				int idx_size = idx.size();
				for(int i = 0; i < idx_size; i++){
					assign_value(str, offsets, idx[i], i);
				}
				n++;
			}
		}
		
		void freeze(){
			if( is_transposed ){
				mat->conservativeResize(Eigen::NoChange, n);
			}else{
				mat->conservativeResize(n, Eigen::NoChange);
			}
		}
		
	private:
		Eigen::MatrixXd* mat;
		bool is_transposed, unknown_dim;
		std::vector<int> idx;
		int n, start_at, chunk_size;
		
		void setSize(int n_fields){
			if( unknown_dim & n == 0){
				int n_cols = n_fields - start_at;
				if( is_transposed ){
					mat->resize(n_cols, n + chunk_size);
				}else{
					mat->resize(n + chunk_size, n_cols);
				}
				for( int i = start_at; i < start_at + n_cols; i++){
					idx.push_back(i);
				}
				unknown_dim = false;
			} 
		}
		
		void checkResize(){
			if( is_transposed ){
				if( mat->cols() < n + 1 ){
					mat->conservativeResize(Eigen::NoChange, n + chunk_size);
				}
			}else{
				if( mat->rows() < n + 1 ){
					mat->conservativeResize(n + chunk_size, Eigen::NoChange);
				}
			}
		}
		
		void assign_value(kstring_t& str, int*& offsets, const int& off_i, const int& mat_i){
			if( is_transposed ){
				mat->coeffRef(mat_i, n) = atof(str.s + offsets[off_i]);
			}else{
				mat->coeffRef(n, mat_i) = atof(str.s + offsets[off_i]);
			}
		}
};

class data_parser
{
	public:
		// Matrix of unknown dimension 
		void add_matrix(Eigen::MatrixXd& new_mat, bool is_trans, int start){
			int nm = v_matrix.size();
			v_matrix.push_back(matrix_filler());
			v_matrix[nm].set(&new_mat, is_trans, start);
		}
		// Matrix of known dimension; keep range of columns 
		void add_matrix(Eigen::MatrixXd& new_mat, const bool is_trans, const int& start, const int& n_fields){
			int nm = v_matrix.size();
			v_matrix.push_back(matrix_filler());
			v_matrix[nm].set(&new_mat, is_trans, start, n_fields);
		}
		// Matrix of known dimension; keep fixed column indices
		void add_matrix(Eigen::MatrixXd& new_mat, const bool is_trans, const std::vector<int>& kp_index, int i_offset = 0){
			int nm = v_matrix.size();
			v_matrix.push_back(matrix_filler());
			v_matrix[nm].set(&new_mat, is_trans, kp_index, i_offset);
		}
		void add_field(std::vector<std::string>& v, const int& i){
			v_string.push_back(&v);
			i_string.push_back(i);
		}
		void add_field(std::vector<long>& v, const int& i){
			v_long.push_back(&v);
			i_long.push_back(i);
		}
		void add_field(std::vector<int>& v, const int& i){
			v_int.push_back(&v);
			i_int.push_back(i);
		}
		void add_field(std::vector<double>& v, const int& i){
			v_double.push_back(&v);
			i_double.push_back(i);
		}
		void add_header(std::vector<std::string>& v, const int& i){
			v_header.push_back(&v);
			i_header.push_back(i);
		}
		
		void add_target(std::vector<std::string>& tg, const int& i, const bool& target_mode_include = true){
			target_keys.push_back(&tg);
			target_cols.push_back(i);
			target_mode.push_back(target_mode_include);
		}
		
		void parse_header(kstring_t& str){
			int n_fields;
			int *offsets = ksplit(&str, 0, &n_fields);
			for( int i = i_header[0]; i < n_fields; i++){
				v_header[0]->push_back(std::string(str.s + offsets[i]));
			}
		}
		
		void parse_fields(kstring_t& str, int*& offsets, int& n_fields)
		{
			if( target_cols.size() > 0 ){
				for( int i = 0; i < target_cols.size(); i++){
					std::string key = std::string(str.s + offsets[target_cols[i]]);
					if( target_mode[i] ){
						if( !has_element(*target_keys[i], key) ){
							return;
						}					
					}else{
						if( has_element(*target_keys[i], key) ){
							return;
						}
					}
				}
			}
			for( int i = 0; i < v_string.size(); i++){
				v_string[i]->push_back(std::string(str.s + offsets[i_string[i]]));
			}
			for( int i = 0; i < v_long.size(); i++){
				v_long[i]->push_back(atol(str.s + offsets[i_long[i]]));
			}
			for( int i = 0; i < v_int.size(); i++){
				v_int[i]->push_back(atoi(str.s + offsets[i_int[i]]));
			}
			for( int i = 0; i < v_double.size(); i++){
				v_double[i]->push_back(atof(str.s + offsets[i_double[i]]));
			}
			for( matrix_filler& mf : v_matrix ){
				mf.push_line(str, offsets, n_fields);
			}
		}
		
		void parse_fields(kstring_t& str)
		{
			int n_fields;
			int *offsets = ksplit(&str, 0, &n_fields);
			parse_fields(str, offsets, n_fields);
		}
		
		template <class hts_file>
		void parse_file(hts_file& htsf, int& n_rows = int0){
			kstring_t str = {0,0,0};
			bool skip_header = ( v_header.size() == 0 );
			while( htsf.next_line(str) >= 0 ){
				if( !str.l ) break;
				if ( str.s[0] == '#' ){
					if( str.s[1] == '#' || (skip_header && n_rows == 0 ) ){
						continue;
					}
				}
				if( n_rows == 0 && i_header.size() > 0 ){
					parse_header(str);
				}else{
					parse_fields(str);
				}
				n_rows++;
			}
			ks_free(&str);
			clear();
		}
		void parse_file(const std::string& fn, int& n_rows = int0)
		{
			basic_hts_file htsf(fn);
			parse_file(htsf, n_rows);
			htsf.close();
		}
		void parse_file(const std::string& fn, const std::string& region, int& n_rows = int0)
		{
			if( region != "" ){
				indexed_hts_file htsf(fn, region);
				parse_file(htsf, n_rows);
				htsf.close();
			}else{
				basic_hts_file htsf(fn);
				parse_file(htsf, n_rows);
				htsf.close();
			}
		}

		void clear(){
			for( matrix_filler& mf : v_matrix ){mf.freeze();}
			v_header.clear();
			i_header.clear();
			v_matrix.clear();
			v_string.clear();
			i_string.clear();
			v_long.clear();
			i_long.clear();
			v_int.clear();
			i_int.clear();
			v_double.clear();
			i_double.clear();
			target_cols.clear();
			target_keys.clear();
			target_mode.clear();
		};

	private:
		std::vector<std::vector<std::string>*> v_header; 
		std::vector<int> i_header;
		std::vector<matrix_filler> v_matrix;
		std::vector<std::vector<std::string>*> v_string;
		std::vector<int> i_string;
		std::vector<std::vector<long>*> v_long;
		std::vector<int> i_long;
		std::vector<std::vector<int>*> v_int;
		std::vector<int> i_int;
		std::vector<std::vector<double>*> v_double;
		std::vector<int> i_double;
		
		std::vector<int> target_cols;
		std::vector<std::vector<std::string>*> target_keys;
		std::vector<bool> target_mode;
};


class file_list
{
	public:
		std::vector<std::string> file_paths;
		
		file_list( const std::string& in_string ){
			
			// Find files from format "path/chr{1:22}.txt"
			std::size_t s_brack = in_string.find("{");
			if( s_brack != std::string::npos ){
				std::string prefix = in_string.substr(0, s_brack);
				std::string suffix = in_string.substr(in_string.find("}")+1);
				
				std::string rstr = in_string.substr(s_brack+1);
				
				int start = std::stoi(rstr.substr(0, rstr.find(":")));
				rstr = rstr.substr(rstr.find(":")+1);
				
				int end = std::stoi(rstr.substr(0, rstr.find("}")));
				
				for(int i = start; i <= end; i++){
					file_paths.push_back(prefix + std::to_string(i) + suffix);
				}
			}else{

				std::vector<std::string> types = {".vcf", ".vcf.gz", ".bcf"};
				
				for( const std::string& tp : types ){
					if( in_string.substr(in_string.size() - tp.size()) == tp ){
						file_paths.push_back(in_string);
					}
				}
					
				if( file_paths.size() == 0 ){
					// No match. Assume that in_string is a list of files.
					data_parser dp;
					dp.add_field(file_paths, 0);
					dp.parse_file(in_string);
				}
			
			}
		};
		
};


#endif