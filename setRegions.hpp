/*
	Corbin's notes on setRegions :
	
	Intersecting variants and regions. May rework logic later. 
	
*/

#ifndef SETREGIONS_HPP
#define SETREGIONS_HPP

#include <vector>

#include "miscUtils.hpp"
#include "setOptions.hpp"
#include "readBED.hpp"
#include "genotypeData.hpp"

using namespace std;

class block_intervals
{
	public:
		vector<int> bed_s;
		vector<int> bed_e;
		vector<int> bcf_s;
		vector<int> bcf_e;
		
		vector<int> bcf_id_s;
		vector<int> bcf_id_e;
		
		vector<int> bed_id_s;
		vector<int> bed_id_e;
		
		vector<string> bcf_regions;
		//vector<string> gene_bcf_regions;
	
		block_intervals(){};	
		block_intervals(bed_data& bdat, genotype_data& gdat, const int& ws){
			make_blocks(bdat, gdat, ws);
		};
		int size(){ return bed_s.size();};
		
		void make_blocks(bed_data&, genotype_data&, const int&);
		
};


#endif




