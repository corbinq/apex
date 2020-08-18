/*
	Corbin's notes on dataParser :
	
	This header provides a generic interface for parsing multiple,  
	heterogeneous fields from htsFiles into multiple, heterogeneous 
	data objects and formats (vector<T> and Eigen::Matrix). 
	
	See readBed.cpp for an example. 
	
*/

#ifndef DATAPARSER_HPP
#define DATAPARSER_HPP

#include <vector>

#include "htsWrappers.hpp"
#include "setOptions.hpp"

using namespace std;

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
				cerr << "Fatal: matrix_filler error, mismatched dimensions\n";
				cerr << n_fields << " != (" << mat->cols() << ", " << mat->rows() << ").\n";
				abort();
			}
		}
		
		// Known matrix size; assign columns from index vector 
		void set(Eigen::MatrixXd* new_mat, const bool is_trans, const vector<int>& kp_index, int i_offset = 0){
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
				cerr << "Fatal: matrix_filler error, mismatched dimensions\n";
				cerr << kp_index.size() << " != (" << mat->cols() << ", " << mat->rows() << ").\n";
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
				int j = 0;
				for(const int& i : idx){
					assign_value(str, offsets, i, j);
					j++;
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
		vector<int> idx;
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
		void add_matrix(Eigen::MatrixXd& new_mat, const bool is_trans, const vector<int>& kp_index, int i_offset = 0){
			int nm = v_matrix.size();
			v_matrix.push_back(matrix_filler());
			v_matrix[nm].set(&new_mat, is_trans, kp_index, i_offset);
		}
		void add_field(vector<string>& v, const int& i){
			v_string.push_back(&v);
			i_string.push_back(i);
		}
		void add_field(vector<long>& v, const int& i){
			v_long.push_back(&v);
			i_long.push_back(i);
		}
		void add_field(vector<int>& v, const int& i){
			v_int.push_back(&v);
			i_int.push_back(i);
		}
		void add_field(vector<double>& v, const int& i){
			v_double.push_back(&v);
			i_double.push_back(i);
		}
		
		void add_target(vector<string>& tg, const int& i, const bool& target_mode_include = true){
			target_keys.push_back(&tg);
			target_cols.push_back(i);
			target_mode.push_back(target_mode_include);
		}
		
		void parse_fields(kstring_t& str, int*& offsets, int& n_fields)
		{
			if( target_cols.size() > 0 ){
				for( int i = 0; i < target_cols.size(); i++){
					string key = string(str.s + offsets[target_cols[i]]);
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
				v_string[i]->push_back(string(str.s + offsets[i_string[i]]));
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
			while( htsf.next_line(str) >= 0 ){
				if( !str.l ) break;
				if ( str.s[0] == '#' && str.s[1] == '#' ) continue;
				parse_fields(str);
				n_rows++;
			}
			ks_free(&str);
			clear();
		}
		void parse_file(const string& fn, int& n_rows = int0)
		{
			basic_hts_file htsf(fn);
			parse_file(htsf, n_rows);
			htsf.close();
		}
		void parse_file(const string& fn, const string& region, int& n_rows = int0)
		{
			indexed_hts_file htsf(fn, region);
			parse_file(htsf, n_rows);
			htsf.close();
		}

		void clear(){
			for( matrix_filler& mf : v_matrix ){mf.freeze();}
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
		vector<matrix_filler> v_matrix;
		vector<vector<string>*> v_string;
		vector<int> i_string;
		vector<vector<long>*> v_long;
		vector<int> i_long;
		vector<vector<int>*> v_int;
		vector<int> i_int;
		vector<vector<double>*> v_double;
		vector<int> i_double;
		
		vector<int> target_cols;
		vector<vector<string>*> target_keys;
		vector<bool> target_mode;
};

#endif