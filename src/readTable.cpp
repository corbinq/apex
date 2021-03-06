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

#include "readTable.hpp"

void table::setRows(std::vector<std::string>& keep){
	rows.setKeepIDs(keep);
	return;
}

void table::setCols(std::vector<std::string>& keep){
	cols.setKeepIDs(keep);
	return;
}

void table::readHeader(const char *file_name)
{
	htsFile *htsf = hts_open(file_name, "r");
	kstring_t str = {0,0,0};
	int line_no = 0;

	while( hts_getline(htsf, KS_SEP_LINE, &str) )
	{
		if( !str.l ) break;
		line_no++;
		if ( str.s[0] == '#' && str.s[1] == '#' ) continue;
		
		int n_fields;
		header_line = line_no;
		int *offsets = ksplit(&str, 0, &n_fields);
		
		// Note: The first column is assumed to be "ID"
		
		for( int i = 1; i < n_fields; i++)
		{
			cols.file.push_back(std::string(str.s + offsets[i]));
		}
		break;
	}
	ks_free(&str);
	hts_close(htsf);
}

void table::readFile(const char *file_name)
{
	if( cols.file.size() == 0 ){
		readHeader(file_name);
	}
	if( cols.keep.size() == 0 ){
		std::vector<std::string> all_cols = cols.file;
		all_cols.erase (all_cols.begin());
		cols.setKeepIDs(all_cols);
	}else if( cols.idx.size() == 0 ){
		cols.makeIndex();
	}
	
	n_cols = cols.keep.size();
	n_rows = 0;
	
	data_matrix = Eigen::MatrixXd(n_cols, 50000);
	
	data_parser dp;
	dp.add_matrix(data_matrix, true, cols.idx, 1);
	
	htsFile *htsf = hts_open(file_name, "r");
	
	kstring_t str = {0,0,0};

	int line_no = 0;

	while( hts_getline(htsf, KS_SEP_LINE, &str) >= 0 )
	{
		if( !str.l ) break;
		
		line_no++;
		if ( str.s[0] == '#' || line_no <= header_line ){
			continue;
		}
		
		int n_fields;
		int* offsets = ksplit(&str, 0, &n_fields);
		
		std::string id = std::string(str.s + offsets[0]);
		
		if ( rows.tryKeep(id) ){
			dp.parse_fields(str, offsets, n_fields);
			n_rows++;
		}
	}
	
	rows.keep = rows.file;

	data_matrix.conservativeResize(Eigen::NoChange, n_rows);
	
	ks_free(&str);
	hts_close(htsf);

	// rows.makeIndex();
}

