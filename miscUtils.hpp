#ifndef MISCUTILS_HPP
#define MISCUTILS_HPP

#include <cmath>
#include <sstream>
#include <vector>
#include <algorithm>

#include "setOptions.hpp"

using namespace std;

static const int print_every_default = 1000;
static vector<int> null_vec = vector<int>(0);

int i_chrom(const string&);
bool ambiguous_snv(const string&, const string&);
string flip_nucleotide(const string&);

vector<string> sort_chroms(vector<string>);

vector<int> seq_int(const int& n);

template<typename T>
inline bool has_element(const vector<T>& v, const T& x){
	return find(v.begin(), v.end(), x) != v.end();
};

void remove_gene_version_number( vector<string>& );

string clean_chrom(const string&);
vector<string> split_string(const string&, const char);

inline void restore_cursor(void){ cerr << "\e[?25h"; };
inline void hide_cursor(void){ cerr << "\e[?25l"; };
inline void clear_line_cerr(void){ cerr << "\33[2K\r"; };

inline void move_back_cerr(const int n){ cerr << string(n, '\b'); };

bool all_lt( const vector<int>& ii, const vector<int>& nn );
bool any_lt( const vector<int>& ii, const vector<int>& nn );
vector<int> which_lt(const vector<int>& ii, const vector<int>& nn);

void print_iter_cerr(int, int, string&);
void thinned_iter_cerr(int& i_last, const int& i_curr, string& suffix, const int& print_every = print_every_default);

void print_header(const vector<string>& cn, ostream& os);

class lindex 
{
	public:
		void set(vector<int>& a){ vals = a; build(); };
		
		lindex() : vals(null_vec), n(-1), a(-1), b(-1) {};
		lindex(vector<int>& a) : vals(a) { build(); };
		
		int index(const int&, bool left = true);
	private:
		int n;
		double a;
		double b;
		vector<int>& vals;
		void build();
};

/*
class block_lindex 
{
	public:
		//lindex();
		block_lindex(vector<string>& chr, vector<int>& pos) : blks(chr),vals(pos) { build(); };
		int index(const string&, const int&);
	private:
		vector<string>& blks;
		vector<int>& vals;
		void build();
};
*/


#endif

