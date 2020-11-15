/*  
    Copyright (C) 2020 
    Author: Corbin Quick <qcorbin@hsph.harvard.edu>

    This file is a part of YAX.

    YAX is distributed "AS IS" in the hope that it will be 
    useful, but WITHOUT ANY WARRANTY; without even the implied 
    warranty of MERCHANTABILITY, NON-INFRINGEMENT, or FITNESS 
    FOR A PARTICULAR PURPOSE.

    The above copyright notice and disclaimer of warranty must 
    be included in all copies or substantial portions of YAX.
*/


/*
	readTable source files and the 'table' class are for processing 
	covariate files, subsetting and filtering covariate columns, and
	subsetting and merging individual IDs in the covariate file with 
	those in the expression file and genotype file. 
	
	HTSLIB's htsFile class is used to seamlessly read from raw text, 
	gzip-compressed, or bgzip compressed input. 
	
	Eigen's MatrixXd class is used to store covariate data. 
	
	Non-numeric discrete covariates (factors) are not currently allowed.
	For now, these should be dummy-coded in the covariate file. 
	
*/

#ifndef READTABLE_HPP
#define READTABLE_HPP

#include <vector>

#include "setOptions.hpp"
#include "dataParser.hpp"
#include "mapID.hpp"

#include <Eigen/Dense>

class table
{
	public:
		int header_line;
	
		int id_column;
		std::string id_column_name;
		
		int n_rows;
		int n_cols;

		id_map rows;
		id_map cols;
		
		Eigen::MatrixXd data_matrix;
		
		void setRows(std::vector<std::string>&);
		void setCols(std::vector<std::string>&);
		void readFile(const char*);
		void readHeader(const char*);
		
		void write_table(const std::string& out_file_path, const std::string& file_header = nullstr){
			
			BGZF* out_file = bgzf_open(out_file_path.c_str(), "w");
			std::stringstream out_line;
			
			
			int dj = 0;
			if( data_matrix.cols() > rows.keep.size() ){
				// has intercept column.
				dj = 1;
			}
			
			
			for(int j = 0; j < rows.keep.size(); j++){
				std::stringstream out_line;
				
				if( j == 0 ){
					if( file_header != nullstr ){
						out_line << file_header << "\n";
					}
					
					out_line << "#id";
					for( const std::string& id : cols.keep ){
						out_line << "\t" << id;
					}
					out_line << "\n";
				}
				
				out_line << rows.keep[j];
				
				for(int i = 0; i < cols.keep.size(); i++){
					out_line << "\t" << data_matrix(i, j + dj);
				}
				out_line << "\n";
				
				write_to_bgzf(out_line.str().c_str(), out_file);
			}
			
			bgzf_close(out_file);
		};
};


#endif

