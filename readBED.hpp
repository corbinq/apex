/*
	Corbin's notes on readBED :
	
	Parsing and storing expression measurements
	in BED file format.
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

using namespace std;

class bed_data
{
	public:
		int header_line;
		
		int n_genes;
	
		id_map ids;
		
		Eigen::MatrixXd data_matrix;
	
		vector<string> chr;
		vector<int> start;
		vector<int> end;
		vector<string> gene_id;
		
		vector<int> block_s;
		vector<int> block_e;
		
		vector<int> v_s;
		vector<int> v_e;
		
		vector<double> stdev;
		
		vector<double> n_var;
		vector<double> pval;
		
		void setKeepIDs(vector<string>& kp_ids){ids.setKeepIDs(kp_ids);};
		
		void readBedFile(const char*, const vector<string>&);
		void readBedHeader(const char*);
};

#endif

