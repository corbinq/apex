# 1 "yax/src/genotypeData.hpp.c"
#ifndef GENOTYPEDATA_HPP
#define GENOTYPEDATA_HPP 

#include <vector>
#include <unordered_map>

#include "setOptions.hpp"
#include "mapID.hpp"
#include "htsWrappers.hpp"
#include "dataParser.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>


void read_sparse_GRM(const std::string& filename, Eigen::SparseMatrix<double>& GRM, const std::vector<std::string>& kp_ids, const double& r_scale, const int& r_col, std::vector<int>& related);

void read_dense_GRM(const std::string&, Eigen::MatrixXd&, std::vector<std::string>&);
# 32 "yax/src/genotypeData.hpp.c"
class sparse_gt
{
 public:
  sparse_gt() : last(-1) {};
  std::vector<std::pair<int, int>> entries;
  int last;
  int size(){
   return entries.size();
  };
  void set_gt(const int& n, const int& gt){
   if( n > last ){
    if( gt != 0 ) entries.push_back(std::make_pair(n,gt));
    last = n;
   }else{
    std::cerr << "Fatal: Tried to add gt out of order!\n";
    exit(2);
   }
  };
  void flip(const int& n){
   int j = 0, i = 0;
   std::vector<std::pair<int, int>> flipped;
   while(i < n){
    if( j < entries.size() ){
     while( i < entries[j].first && i < n ){
      flipped.push_back( std::make_pair(i, 2) );
      i++;
     }
     if( i == entries[j].first ){
      if( entries[j].second == 1 ){
       flipped.push_back(entries[j]);
      }else if(entries[j].second != 2){
       std::cerr <<"\n\n" << entries[j].second << "\n\n";
       std::cerr << "entries[j].second != 2\n";
       exit(1);
      }
      i++;
      j++;
     }else if( j > i ){
      std::cerr << "This should never happen!\n";
      exit(2);
     }
    }else{
     flipped.push_back( std::make_pair(i, 2) );
     i++;
    }
   }
   entries = flipped;
  };
  void add_gt_sparsemat(Eigen::SparseMatrix<double>& smat, const int& col_n){
   smat.startVec(col_n);
   for( const auto& x : entries){
    smat.insertBack(x.first, col_n) = x.second;
   }
  };
};


class sparse_ds
{
 public:
  std::vector<float> dosages;
  int last;

  sparse_ds() : last(-1) {};
  sparse_ds(const int& N) : last(-1) {dosages.resize(N);};

  void set_ds(const int& n, const float& ds){
   if( n > last ){
    dosages[n] = ds;
    last = n;
   }else{
    std::cerr << "Fatal: Tried to add ds out of order!\n";
    exit(2);
   }
  };
  void flip(){
   for( float& ds: dosages ){
    ds = 2 - ds;
   }
  };
  void add_ds_sparsemat(Eigen::SparseMatrix<double>& smat, const int& col_n){
   smat.startVec(col_n);
   for(int i = 0; i < dosages.size(); i++){
    if( dosages[i] > 0 ){
     smat.insertBack(i, col_n) = (double) dosages[i];
    }
   }
  };
};


class genotype_data
{
 public:
  genotype_data(): chunk_size(50000), nz_frac(0.08), n_variants(0), n_samples(0) {};

  int chunk_size;
  double nz_frac;

  int n_samples;
  int n_variants;

  id_map ids;

  std::vector<std::string> chr;
  std::vector<int> pos;

  std::vector<std::string> ref;
  std::vector<std::string> alt;

  std::vector<std::pair<int,int>> ld_index_range;

  std::vector<double> mean;
  std::vector<double> var;

  std::vector<bool> flipped;
  std::vector<bool> keep;

  int geno_start;
  int geno_size;

  Eigen::SparseMatrix<double> genotypes;

  void read_bcf_variants(bcf_srs_t*&, bcf_hdr_t*&, int&, bool store_geno = true, bool scan_geno = true);

  int get_ld_index();
  inline bool process_bcf_variant(bcf1_t*&, bcf_hdr_t*&, bool store_geno = true, bool scan_geno = true);
  void read_bcf_header(bcf_hdr_t*);
  void freeze_genotypes();
  void resize_genotypes(const int& nv){
   genotypes.resize(ids.keep.size(), nv);
   genotypes.reserve( (int) (nz_frac * genotypes.rows() * genotypes.cols() ) );
  };
  void initialize_genotypes(const int& nv){
   resize_genotypes(nv);
  };
  void print_genotypes(){
   std::cout << genotypes << "\n";
  };

  void clear();
  void clear_genotypes();
  void read_genotypes(bcf_srs_t*&, bcf_hdr_t*&, const int&, const int&);

  bool record_matches(bcf1_t*& rec, bcf_hdr_t*& hdr, const int& i);

 private:
  inline bool add_bcf_genotypes(int*&, const int&, double&, double&, bool&, const bool);
  inline bool add_bcf_dosages(float*& ds_rec, const int& col_n, double& mean_, double& var_, bool& flip_, const bool store_geno);

  void checkResize(const int& col_n){
   if( col_n >= genotypes.cols() ){
    genotypes.conservativeResize(ids.keep.size(), col_n + chunk_size);
   }
  };

};


#endif
