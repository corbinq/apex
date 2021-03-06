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

/*
mapID source files and the 'id_map' class are for subsetting, merging, 
and mapping row and column IDs in across multiple files (expression, 
covariate, and genotype). These classes help reduce redundancy for 
common tasks across file types. 	
*/

#ifndef MAPID_HPP
#define MAPID_HPP

#include <vector>
#include <algorithm>
#include <unordered_set>

#include "miscUtils.hpp"
#include "setOptions.hpp"


std::vector<std::string> intersect_ids(std::vector<std::string>, std::vector<std::string>);

class id_map
{
	public:
		std::unordered_set<std::string> keep_set;
	
		std::vector<std::string> file;
		std::vector<std::string> keep;
		std::vector<int> idx;
		std::vector<int> idx_f2k;
		
		void makeIndex();
		bool tryKeep(std::string&);
		void setFileIDs(std::vector<std::string>&);
		void setKeepIDs(std::vector<std::string>&);
		
		int n();
};

#endif




