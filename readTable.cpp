/*

	Corbin's notes on readTable :
	
	COVARIATE DATA
	
	readTable source files and the 'table' class are for processing 
	covariate files, subsetting and filtering covariate columns, and
	subsetting and merging individual IDs in the covariate file with 
	those in the expression file and genotype file. 
	
	We need to take some care in the file format. Some covariate files
	store each sample as a row, while others store each covariate as 
	a row (similar to expression bed files)
	
	HTSLIB's htsFile class is used to seamlessly read from raw text, 
	gzip-compressed, or bgzip compressed input. 
	
	Eigen's MatrixXd class is used to store covariate data. 
	
	Non-numeric discrete covariates (factors) are not currently allowed.
	For now, these should be dummy-coded in the covariate file. 
	
*/

#include "readTable.hpp"

using namespace std;

void table::setRows(vector<string>& keep){
	rows.setKeepIDs(keep);
	return;
}

void table::setCols(vector<string>& keep){
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
			cols.file.push_back(string(str.s + offsets[i]));
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
		vector<string> all_cols = cols.file;
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
		
		string id = string(str.s + offsets[0]);
		
		if ( rows.tryKeep(id) ){
			dp.parse_fields(str, offsets, n_fields);
			n_rows++;
		}
	}

	data_matrix.conservativeResize(Eigen::NoChange, n_rows);
	
	ks_free(&str);
	hts_close(htsf);

	// rows.makeIndex();
}

