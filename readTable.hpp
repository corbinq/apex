/*
	Corbin's notes on readTable :
	
	COVARIATE DATA
	
	readTable source files and the 'table' class are for processing 
	covariate files, subsetting and filtering covariate columns, and
	subsetting and merging individual IDs in the covariate file with 
	those in the expression file and genotype file. 
	
	HTSLIB's htsFile class is used to seamlessly read from raw text, 
	gzip-compressed, or bgzip compressed input. 
	
	Eigen's MatrixXd class is used to store covariate data. 
	
	Non-numeric discrete covariates (factors) are not currently allowed.
	For now, these should be dummy-coded in the covariate file. 
	
*/

#ifndef READTABLE_HPP
#define READTABLE_HPP

#include <vector>

#include "setOptions.hpp"
#include "dataParser.hpp"
#include "mapID.hpp"

#include <Eigen/Dense>

using namespace std;

class table
{
	public:
		int header_line;
	
		int id_column;
		string id_column_name;
		
		int n_rows;
		int n_cols;

		id_map rows;
		id_map cols;
		
		Eigen::MatrixXd data_matrix;
		
		void setRows(vector<string>&);
		void setCols(vector<string>&);
		void readFile(const char*);
		void readHeader(const char*);
};


#endif

