# 1 "yax/src/processVCOV.hpp.c"
#ifndef PROCESSVCOV_HPP
#define PROCESSVCOV_HPP 

#include <vector>

#include "setOptions.hpp"
#include "readTable.hpp"
#include "genotypeData.hpp"
#include "fitUtils.hpp"
#include "miscUtils.hpp"
#include "dataParser.hpp"
#include "xzFormat.hpp"

#include <shrinkwrap/istream.hpp>
#include <Eigen/Sparse>
#include <Eigen/Dense>


typedef uint16_t vbin_t;
static const long double VBIN_T_MAX = (pow(2,8*sizeof(vbin_t)) - 1);

static const bool xz_mode = true;
static const size_t XZ_BUFFER_SIZE = 2 * 1024 * 1024;

class vcov_bin_file
{
 public:
  std::string file_name;

  BGZF *fp;
  xzReader fp_xz;

  bool is_xz;

  vcov_bin_file() {};
  vcov_bin_file(std::string);
  void open(std::string);
  void close();
  std::vector<vbin_t> get_data(long,long);
};

class vcov_data
{
 public:
  std::string file_prefix;
  std::string region;

  std::vector<std::string> chr;
  std::vector<int> pos;

  std::vector<std::string> ref;
  std::vector<std::string> alt;
  std::vector<int> flipped;

  std::vector<double> mac;
  std::vector<double> var;

  std::vector<int> rid;
  std::vector<long> bin_start_idx;
  std::vector<long> bin_nvals;

  vcov_bin_file bin;

  Eigen::MatrixXd GtU;

  vcov_data() {};
  vcov_data(std::string,std::string);

  void open(std::string,std::string);
  void close();

  Eigen::MatrixXd getGtG(const std::vector<int>& v_i, const std::vector<int>& v_j = std::vector<int>(0));
  Eigen::MatrixXd getV(const std::vector<int>& v_i, const std::vector<int>& v_j = std::vector<int>(0));

  Eigen::MatrixXd getGtG(int, int);
  Eigen::MatrixXd getV(int, int);

  Eigen::MatrixXd getColGtG(int, int, int);
  Eigen::MatrixXd getColV(int, int, int);

  Eigen::MatrixXd getColGtG(int, const std::vector<int>& );
  Eigen::MatrixXd getColV(int, const std::vector<int>& );

  double getPairV(int, int);

  Eigen::MatrixXd getV_i(const std::vector<int>&, const int offset = 0);
};

class vcov_bin_gz_out
{
 public:
  vcov_bin_gz_out(std::string fn) :
   file_name(fn),
   fp_xz(fn.c_str())
  {
   file_name = fn;
   bytes_c = 0;

   if( xz_mode ){
    buf_size = XZ_BUFFER_SIZE;
    std::cerr << "\nWriting data to " << file_name << "\n";
   }else{
    buf_size = BUFFER_SIZE;
    std::string file_name_gz = file_name + ".gz";
    fp = bgzf_open(file_name_gz.c_str(), "w\0");

    bgzf_index_build_init(fp);
   }
  };

  void open(std::string);
  void add_to_queue(vbin_t);
  void close();

 private:
  std::vector<vbin_t> queue;

  size_t buf_size;

  BGZF *fp;
  shrinkwrap::xz::ostream fp_xz;

  std::string file_name;

  void write_queue();

  size_t bytes_c;
};


void read_LD_gz_bytes(std::string);

void write_vcov_files(genotype_data&, const table&);

#endif
