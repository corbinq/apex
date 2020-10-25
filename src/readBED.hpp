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
	
		std::vector<std::string> chr;
		std::vector<int> start;
		std::vector<int> end;
		std::vector<std::string> gene_id;
		
		std::vector<int> block_s;
		std::vector<int> block_e;
		
		std::vector<int> v_s;
		std::vector<int> v_e;
		
		std::vector<double> stdev;
		
		std::vector<double> n_var;
		std::vector<double> pval;
		
		void setKeepIDs(std::vector<std::string>& kp_ids){ids.setKeepIDs(kp_ids);};
		
		void readBedFile(const char*, const std::vector<std::string>&);
		void readBedHeader(const char*);
};

#endif

