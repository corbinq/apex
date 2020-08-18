#include "miscUtils.hpp"

using namespace std;

bool all_lt( const vector<int>& ii, const vector<int>& nn ){
	for(int i = 0; i < ii.size(); ++i){
		if( ii[i] >= nn[i] ){
			return false;
		}
	}
	return true;
}

bool any_lt( const vector<int>& ii, const vector<int>& nn ){
	for(int i = 0; i < ii.size(); ++i){
		if( ii[i] < nn[i] ){
			return true;
		}
	}
	return false;
}

vector<int> which_lt(const vector<int>& ii, const vector<int>& nn){
	vector<int> out;
	for(int i = 0; i < ii.size(); ++i){
		if( ii[i] < nn[i] ){
			out.push_back(i);
		}
	}
	return out;
}


void remove_gene_version_number( vector<string>& ids){
	for( string& id : ids ){
		 id = id.substr(0, id.find("."));
	}
	return;
}

vector<int> seq_int(const int& n){
	vector<int> out;
	for(int i = 0; i < n; i++){
		out.push_back(i);
	}
	return out;
}

vector<string> split_string(const string& input, const char delim)
{
	stringstream ss(input);
	string field;
	vector<string> out;
	
	while( getline(ss, field, delim) ){
		out.push_back(field);
	}
	return out;
}
		
string clean_chrom(const string& x){
	if( x.length() >= 3 ){
		if( x.substr(0,3) == "chr" || x.substr(0,3) == "CHR" ){
			return x.substr(3, x.length() - 3);
		}
	}
	return x;
}

int i_chrom( const string& chrom )
{
	auto f_chrom = global_opts::i_chrom_map.find(chrom);
	if( f_chrom == global_opts::i_chrom_map.end() ){
		return -1;
	}else{
		return f_chrom->second;
	}
}

vector<string> sort_chroms( vector<string> chroms ){
	sort(chroms.begin(), chroms.end(), 
          [](const string& x, const string& y){ return i_chrom(x) < i_chrom(y); }
	);
	return chroms;
}


bool ambiguous_snv( const string& ref, const string& alt )
{
	if( ref == "T" || ref == "G" ){
		return ambiguous_snv(alt, ref);
	}else if( ref == "A" && alt == "T" ){
		return true; 
	}else if( ref == "C" && alt == "G" ){
		return true;
	}
	return false;
}

string flip_nucleotide(const string& ref)
{
	if( ref == "A" ) return "T";
	if( ref == "T" ) return "A";
	if( ref == "C" ) return "G";
	if( ref == "G" ) return "C";
	return "";
}

void print_iter_cerr(int i_last, int i_curr, string& suffix){
	string s_last = to_string(i_last);
	string s_curr = to_string(i_curr);
	if( i_last < i_curr ){
		move_back_cerr(s_last.length());
	}
	cerr << i_curr;
	if( s_last.length() != s_curr.length() || i_last > i_curr ){
		cerr << suffix;
		move_back_cerr(suffix.length());
	}
	return;
}

void thinned_iter_cerr(int& i_last, const int& i_curr, string& suffix, const int& print_every){
	if( rand() % print_every == 0 ){
		string s_last = to_string(i_last);
		string s_curr = to_string(i_curr);
		if( i_last < i_curr ){
			move_back_cerr(s_last.length());
		}
		cerr << i_curr;
		if( s_last.length() != s_curr.length() || i_last > i_curr ){
			cerr << suffix;
			move_back_cerr(suffix.length());
		}
		i_last = i_curr;
	}
	return;
}

void print_header(const vector<string>& cn, ostream& os){
	bool is_first = true;
	for( const string& s : cn ){
		if( is_first ){
			is_first = false;
		}else{
			os << "\t";
		}
		os << s;
	}
	os << "\n";
	return;
}

void lindex::build(){
	
	double m_xy = 0;
	double m_x = 0;
	double m_xx = 0;
	double m_y = 0;
	n = 0;
	
	double sc = 0.000001;
	for( const int& x : vals){
		double xs = x * sc;
		m_x += xs;
		m_xx += xs*xs;
		m_y += n;
		m_xy += n*xs;
		n++;
	}
	m_y /= n;
	m_x /= n;
	m_xy /= n;
	m_xx /= n;
	
	b = ( m_xy - m_x*m_y )/( m_xx - m_x*m_x );
	a = m_y - b*m_x;
	b = b * sc;
}

int lindex::index(const int& x, bool left){
	int i = round(a + x*b);
	int first = i;
	i = i < n ? i : n-1;
	i = i > 0 ? i : 0;
	
	if( left ){
		while( vals[i] >= x && i > 0 ){
			i--;
		}
		while( vals[i] < x ){
			if( i + 1 < n ){
				if( vals[i + 1] <= x ){
					i++;
				}else{
					break;
				}
			}else{
				break;
			}
		}
	}else{
		while( vals[i] > x && i > 0 ){
			i--;
		}
		while( vals[i] <= x ){
			if( i + 1 < n ){
				if( vals[i + 1] <= x ){
					i++;
				}else{
					break;
				}
			}else{
				break;
			}
		}
	}
	
	return i;
}

