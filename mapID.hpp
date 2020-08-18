/*
	Corbin's notes on mapID :
	
	ROW/COLUMN ID MANAGEMENT
	
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

using namespace std;

vector<string> intersect_ids(vector<string>, vector<string>);

class id_map
{
	public:
		unordered_set<string> keep_set;
	
		vector<string> file;
		vector<string> keep;
		vector<int> idx;
		vector<int> idx_f2k;
		
		void makeIndex();
		bool tryKeep(string&);
		void setFileIDs(vector<string>&);
		void setKeepIDs(vector<string>&);
		
		int n();
};

#endif




