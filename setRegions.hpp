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


class block_intervals
{
	public:
		std::vector<int> bed_s;
		std::vector<int> bed_e;
		std::vector<int> bcf_s;
		std::vector<int> bcf_e;
		
		std::vector<int> bcf_id_s;
		std::vector<int> bcf_id_e;
		
		std::vector<int> bed_id_s;
		std::vector<int> bed_id_e;
		
		std::vector<std::string> bcf_regions;
		//std::vector<std::string> gene_bcf_regions;
	
		block_intervals(){};	
		block_intervals(bed_data& bdat, genotype_data& gdat, const int& ws){
			make_blocks(bdat, gdat, ws);
		};
		int size(){ return bed_s.size();};
		
		void make_blocks(bed_data&, genotype_data&, const int&);
		
};


#endif




