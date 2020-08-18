#ifndef PROCESSVCOV_HPP
#define PROCESSVCOV_HPP

#include <vector>

#include "setOptions.hpp"
#include "readTable.hpp"
#include "genotypeData.hpp"
#include "fitUtils.hpp"
#include "miscUtils.hpp"
#include "dataParser.hpp"

#include <Eigen/Sparse>
#include <Eigen/Dense>

using namespace std;

typedef uint16_t vbin_t;
static const long double VBIN_T_MAX = (pow(2,8*sizeof(vbin_t)) - 1);

class vcov_bin_file
{
	public:
		string file_name;
		
		BGZF *fp;
		
		vcov_bin_file() {};
		vcov_bin_file(string);
		void open(string);
		void close();
		vector<vbin_t> get_data(long,long);
};

class vcov_data
{
	public:
		string file_prefix;
		string region;
	
		vector<string> chr;
		vector<int> pos;
		//vector<string> rsid; 
		vector<string> ref;
		vector<string> alt;
		vector<int> flipped;
		
		vector<double> mac;
		vector<double> var;

		vector<int> rid;
		vector<long> bin_start_idx;
		vector<long> bin_nvals;
		
		vcov_bin_file bin;
		
		Eigen::MatrixXd GtU;

		vcov_data() {};
		vcov_data(string,string);
		
		void open(string,string);
		void close();
		
		Eigen::MatrixXd getGtG(const vector<int>& v_i, const vector<int>& v_j = vector<int>(0));
		Eigen::MatrixXd getV(const vector<int>& v_i, const vector<int>& v_j = vector<int>(0));
		
		Eigen::MatrixXd getGtG(int, int);
		Eigen::MatrixXd getV(int, int);
		
		Eigen::MatrixXd getColGtG(int, int, int);
		Eigen::MatrixXd getColV(int, int, int);
		
		Eigen::MatrixXd getColGtG(int, const vector<int>& );
		Eigen::MatrixXd getColV(int, const vector<int>& );
		
		double getPairV(int, int);
		
		Eigen::MatrixXd getV_i(const vector<int>&, const int offset = 0);
};

class vcov_bin_gz_out
{
	public:
		vcov_bin_gz_out(string);
		
		void open(string);
		void add_to_queue(vbin_t);
		void close();
		
	private:
		vector<vbin_t> queue;
		
		BGZF *fp;
		void *buffer;
		
		string file_name;

		void write_queue();
};


void read_LD_gz_bytes(string);

void write_vcov_files(genotype_data&, const table&);

#endif

