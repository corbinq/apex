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




