/*  
    Copyright (C) 2020 
    Author: Corbin Quick <qcorbin@hsph.harvard.edu>

    This file is a part of APEX.

    APEX is distributed "AS IS" in the hope that it will be 
    useful, but WITHOUT ANY WARRANTY; without even the implied 
    warranty of MERCHANTABILITY, NON-INFRINGEMENT, or FITNESS 
    FOR A PARTICULAR PURPOSE.

    The above copyright notice and disclaimer of warranty must 
    be included in all copies or substantial portions of APEX.
*/

#ifndef READBED_HPP
#define READBED_HPP

#include <iostream>
#include <stdexcept>
#include <vector>
#include <unordered_map>
#include <algorithm>

#include "setOptions.hpp"
#include "mapID.hpp"
#include "dataParser.hpp"

#include <Eigen/Dense>

class bed_data
{
	public:
		int header_line;
		
		int n_genes;
	
		id_map ids;
		
		Eigen::MatrixXd data_matrix;
		
		bool is_residual;
	
		std::vector<std::string> chr;
		std::vector<int> start;
		std::vector<int> end;
		std::vector<std::string> gene_id;
		
		std::vector<int> block_s;
		std::vector<int> block_e;
		
		std::vector<int> v_s;
		std::vector<int> v_e;

    // Vector of standard deviations for each gene
		std::vector<double> stdev;
		
		std::vector<double> n_var;
		std::vector<double> pval;
		
		void setKeepIDs(std::vector<std::string>& kp_ids){ids.setKeepIDs(kp_ids);};
		
		void readBedFile(const char*, const std::vector<std::string>&);
		void readBedHeader(const char*);
		
		void write_bed(const std::string& out_file_path, const std::string& bed_header = nullstr){
			
			BGZF* out_file = bgzf_open(out_file_path.c_str(), "w");
			std::stringstream out_line;
			
			for(int j = 0; j < gene_id.size(); j++){
				std::stringstream out_line;
				
				if( j == 0 ){
					if( bed_header != nullstr ){
						out_line << bed_header << "\n";
					}
					
					out_line << "#chr" << "\t" << "start" << "\t" << "end" << "\t" << "gene";
					for( const std::string& id : ids.keep ){
						out_line << "\t" << id;
					}
					out_line << "\n";
				}
				
				out_line << chr[j] << "\t" << start[j] << "\t" << end[j] << "\t" << gene_id[j];
				
				for(int i = 0; i < ids.keep.size(); i++){
					out_line << "\t" << data_matrix(i,j);
				}
				out_line << "\n";
				
				write_to_bgzf(out_line.str().c_str(), out_file);
			}
			
			bgzf_close(out_file);
			build_tabix_index(out_file_path, 1);
			
		};
};

#endif

