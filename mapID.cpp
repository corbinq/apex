/*
	Corbin's notes on mapID :
	
	ROW/COLUMN ID MANAGEMENT
	
	mapID source files and the 'id_map' class are for subsetting, merging, 
	and mapping row and column IDs in across multiple files (expression, 
	covariate, and genotype). These classes help reduce redundancy for 
	common tasks across file types. 
	
*/

#include "mapID.hpp"

using namespace std;

vector<string> intersect_ids(vector<string> a, vector<string> b)
{
	vector<string> out;
	sort(a.begin(), a.end());
	sort(b.begin(), b.end());
	set_intersection(a.begin(), a.end(), b.begin(), b.end(), back_inserter(out));
	return out;
}

int id_map::n()
{	
	return keep.size();
}

void id_map::setFileIDs(vector<string>& ids)
{	
	file = ids;
	
	if( keep.size() > 0 ){
		makeIndex();
	}
}

bool id_map::tryKeep(string& id)
{
	bool should_keep = true;
	if( keep.size() > 0 ){
		if( keep_set.size() == 0 ){
			copy(keep.begin(),keep.end(),inserter(keep_set,keep_set.end()));
		}
		should_keep = keep_set.find(id) != keep_set.end();
	}
	if( should_keep ){
		file.push_back(id);
	}
	return should_keep;
}

void id_map::makeIndex()
{
	idx.clear();
	idx_f2k.clear();
	
	unordered_map<string, int> file_id_map;
	unordered_map<string, int> keep_id_map;
	
	int ii = 0;
	for( string& id : file )
	{
		file_id_map[id] = ii;
		ii++;
	}
	ii = 0;
	vector<int> not_in_file;
	for( string& id : keep )
	{
		if( file_id_map.find(id) == file_id_map.end() ){
			not_in_file.push_back(ii);
		}else{
			idx.push_back(file_id_map[id]);
		}
		ii++;
	}
	for(auto j = not_in_file.rbegin(); j != not_in_file.rend(); ++j)
	{ 
		keep.erase(keep.begin() + *j);
	}
	
	
	ii = 0;
	for( string& id : keep )
	{
		keep_id_map[id] = ii;
		ii++;
	}
	for( string& id : file )
	{
		if( keep_id_map.find(id) == keep_id_map.end() ){
			idx_f2k.push_back(-1);
		}else{
			idx_f2k.push_back(keep_id_map[id]);
		}
	}
}

void id_map::setKeepIDs(vector<string>& ids)
{	
	keep = ids;
	
	// map the ids we want to keep to the ids in the file
	if( file.size() > 0 )
	{	
		makeIndex();
	}else
	{
		// if we set 'keep' without setting file ids, we will make
		// a hash to check whether each file id should be kept. 
		copy(keep.begin(),keep.end(),inserter(keep_set,keep_set.end()));
	}
}

