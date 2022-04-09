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


#ifndef MISCUTILS_HPP
#define MISCUTILS_HPP

#include <cmath>
#include <sstream>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <string>
#include <memory>
#include "setOptions.hpp"

template<typename ... Args>
std::string string_format(const std::string& format, Args ... args) {
  int size_s = std::snprintf(nullptr, 0, format.c_str(), args ...) + 1;
  if (size_s <= 0) {
    throw std::runtime_error("Error during formatting.");
  }

  auto size = static_cast<size_t>(size_s);
  auto buf = std::make_unique<char[]>(size);
  std::snprintf(buf.get(), size, format.c_str(), args ...);
  return {buf.get(), buf.get() + size - 1};
}

template <typename float_type> std::string log_to_string(float_type& value) {
  if (isnan(value)) {
    return "nan";
  }
  else if (isinf(value)) {
    if (signbit(value)) {
      return "-inf";
    }
    else {
      return "inf";
    }
  }

  float_type v = value / M_LN10;
  int exp = floor(v);
  float_type r = v - exp;
  float_type base = pow(10, r);
  std::string text = string_format("%fe%i", base, exp);
  return text;
}

static const int print_every_default = 1000;
static std::vector<int> null_vec = std::vector<int>(0);

int i_chrom(const std::string&);
bool ambiguous_snv(const std::string&, const std::string&);
std::string flip_nucleotide(const std::string&);

std::vector<std::string> sort_chroms(std::vector<std::string>);

std::vector<int> seq_int(const int& n);

template<typename T>
inline bool has_element(const std::vector<T>& v, const T& x){
	return find(v.begin(), v.end(), x) != v.end();
};

template<typename T>
inline bool has_element(const std::unordered_set<T>& v, const T& x){
	return v.find(x) != v.end();
};

void remove_gene_version_number( std::vector<std::string>& );

std::string clean_chrom(const std::string&);
std::vector<std::string> split_string(const std::string&, const char);

inline void restore_cursor(void){ std::cerr << "\e[?25h"; };
inline void hide_cursor(void){ std::cerr << "\e[?25l"; };
inline void clear_line_cerr(void){ std::cerr << "\33[2K\r"; };

inline void move_back_cerr(const int n){ std::cerr << std::string(n, '\b'); };

bool all_lt( const std::vector<int>& ii, const std::vector<int>& nn );
bool any_lt( const std::vector<int>& ii, const std::vector<int>& nn );
std::vector<int> which_lt(const std::vector<int>& ii, const std::vector<int>& nn);

void print_iter_cerr(int, int, std::string&);
void thinned_iter_cerr(int& i_last, const int& i_curr, std::string& suffix, const int& print_every = print_every_default);

void print_header(const std::vector<std::string>& cn, std::ostream& os);

class lindex 
{
	public:
		void set(std::vector<int>& a){ vals = a; build(); };
		
		lindex() : vals(null_vec), n(-1), a(-1), b(-1) {};
		lindex(std::vector<int>& a) : vals(a) { build(); };
		
		int index(const int&, bool left = true);
	private:
		int n;
		double a;
		double b;
		std::vector<int>& vals;
		void build();
};

/*
class block_lindex 
{
	public:
		//lindex();
		block_lindex(std::vector<std::string>& chr, std::vector<int>& pos) : blks(chr),vals(pos) { build(); };
		int index(const std::string&, const int&);
	private:
		std::vector<std::string>& blks;
		std::vector<int>& vals;
		void build();
};
*/


#endif

