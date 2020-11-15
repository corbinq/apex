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


/*
Generic wrappers for HTSLIB BCF/VCF file parsing, as well as for BGZF reading/writing.
*/

#ifndef HTSWRAPPERS_HPP
#define HTSWRAPPERS_HPP

#include <vector>
#include <cstdlib>

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "htslib/tbx.h"
#include "htslib/hts.h"
#include "htslib/bgzf.h"
#include "htslib/kseq.h"
#include "htslib/kstring.h"
#include "htslib/synced_bcf_reader.h"

#include <Eigen/Dense>

#include "setOptions.hpp"
#include "miscUtils.hpp"


static const int BUFFER_SIZE = 64 * 1024;

static int MINUS99 = -99;

static const std::string null_string = "";

static Eigen::MatrixXd NULL_MATRIX = Eigen::MatrixXd::Zero(0,0);

static std::vector<int> null_vector_int = std::vector<int>(0);

class basic_hts_file
{
	public:
		void open(const std::string& prefix){htsf = hts_open(prefix.c_str(), "r");};
		int next_line(kstring_t& str){ return hts_getline(htsf, KS_SEP_LINE, &str); };

		basic_hts_file(){};
		basic_hts_file(const std::string& prefix){ open(prefix); };
		
		void close(){ hts_close(htsf);};
		
	private:
		htsFile *htsf;
};

class indexed_hts_file
{
	public:
	
		void open(const std::string& prefix, const std::vector<std::string>& reg);
		void open(const std::string& prefix, const std::string& reg = null_string){
			std::vector<std::string> tmp_reg;
			if( reg != "" ) tmp_reg.push_back(reg);
			open(prefix, tmp_reg);
		};

		int next_line(kstring_t& str);

		indexed_hts_file(){};
		indexed_hts_file(const std::string& prefix, const std::string& reg = null_string){
			open(prefix, reg);
		};
		indexed_hts_file(const std::string& prefix, const std::vector<std::string>& reg){
			open(prefix, reg);
		};

		void close(){
			if (itr) hts_itr_destroy(itr);
			hts_close(htsf);
		};
		
	private:
		int n_reg;
		int c_reg;
		std::vector<std::string> regions;
		
		htsFile *htsf;
		tbx_t *tbx;
		hts_itr_t *itr;
};

void bcf_seek(bcf_srs_t*, const std::string&, const int& start = MINUS99, const int& end = MINUS99);

void bgzip_file(std::string, int);

inline void write_to_bgzf(const std::string& inp, BGZF *fp){
	if( bgzf_write(fp, (const void*) inp.c_str(), inp.length()) < 0 ){
		std::cerr << "Error: Could not write to output file\n";
		exit(1);
	}
};

inline void write_to_bgzf(const char* inp, BGZF *fp){
	if( bgzf_write(fp, (const void*) inp, strlen(inp)) < 0 ){
		std::cerr << "Error: Could not write to output file\n";
		exit(1);
	}
};

inline void build_tabix_index(std::string file_name, int type = 0){
	if(tbx_index_build(file_name.c_str(), 14, type==0? &tbx_conf_vcf : &tbx_conf_bed )!=0)
	{
		std::cerr << "Fatal error: Couldn't create index for " << file_name << "\n";
		exit(1);
	}
};

std::vector<std::string> get_chroms(std::string file_name, std::vector<int>& n_variants = null_vector_int);

#endif

