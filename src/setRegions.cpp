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



#include "setRegions.hpp"


void seek_to(int& i, const std::vector<std::string>& chr, const std::vector<int>& pos, const std::string& target_chr, const int target_pos){
	while( chr[i] != target_chr && i < chr.size() - 1 ){
		i++;
	}
	while( chr[i] == target_chr && pos[i] < target_pos && i < chr.size() - 1 ){
		if( chr[i + 1] == target_chr && pos[i + 1]  < target_pos ){
			i++;
		}else{
			break;
		}
	}
}

void block_intervals::make_blocks(bed_data& bdat, genotype_data& gdat, const int& ws){

	std::vector<std::string>& chr_b = bdat.chr;
	std::vector<int>& start_b = bdat.start; 
	std::vector<int>& end_b = bdat.end; 
	
	std::vector<std::string>& chr_g = gdat.chr; 
	std::vector<int>& pos_g = gdat.pos;

	std::vector<int>& b_s_b = bdat.block_s;
	std::vector<int>& b_e_b = bdat.block_e;
	
	std::vector<int>& v_s_b = bdat.v_s;
	std::vector<int>& v_e_b = bdat.v_e;
	

	bool block_by_variant = false; 

	// int block_size = bs;
	int window = ws;
	
	int block_id = 0;
	
	int n_b = chr_b.size();
	int n_g = chr_g.size();
	
	bcf_id_s = std::vector<int>(n_g,-1);
	bcf_id_e = std::vector<int>(n_g,-1);
	
	bed_id_s = std::vector<int>(n_b,-1);
	bed_id_e = std::vector<int>(n_b,-1);
	
	b_s_b = std::vector<int>(n_b,-1);
	b_e_b = std::vector<int>(n_b,-1);
	
	v_s_b = std::vector<int>(n_b,-1);
	v_e_b = std::vector<int>(n_b,-1);
	
	int i_b = 0; int i_g = 0;
	int i_b_e = 0; int i_g_e = 0;
	
	double magic_fraction = 0.8;
	
	while( i_b < n_b ){
		
		i_b_e = i_b;
		
		seek_to(i_b_e, chr_b, start_b, chr_b[i_b], end_b[i_b] + magic_fraction * window);
		seek_to(i_g, chr_g, pos_g, chr_b[i_b], start_b[i_b] - window);
		
		if( chr_g[i_g] == chr_b[i_b] ){
			
			i_g_e = i_g_e > i_g ? i_g_e : i_g;
			seek_to(i_g_e, chr_g, pos_g, chr_b[i_b], end_b[i_b_e] + window + 1);
			
			if( i_b_e >= i_b && pos_g[i_g] - end_b[i_b] <= window && chr_g[i_g] == chr_b[i_b] ){
				
				bed_s.push_back(i_b);
				bed_e.push_back(i_b_e);
				bcf_s.push_back(i_g);
				bcf_e.push_back(i_g_e);
				
				for( int j = i_g; j <= i_g_e; ++j ){
					if( bcf_id_s[j] < 0 ){
						bcf_id_s[j] = block_id;
					}
					bcf_id_e[j] = block_id;
				}
				
				std::string region = chr_g[i_g] + ":" + std::to_string(pos_g[i_g]) + "-" + std::to_string(pos_g[i_g_e]);
				
				bcf_regions.push_back(region);
	
				for( int i = i_b; i <= i_b_e; ++i ){
					if( bed_id_s[i] < 0 ){
						bed_id_s[i] = block_id;
						b_s_b[i] = block_id;
					}
					bed_id_e[i] = block_id;
					b_e_b[i] = block_id;
					
					v_s_b[i] = (i == i_b ? i_g : v_s_b[i-1]);
					seek_to(v_s_b[i], chr_g, pos_g, chr_b[i_b], start_b[i] - window);
					
					v_e_b[i] = (i == i_b ? i_g : v_e_b[i-1]);
					
					seek_to(v_e_b[i], chr_g, pos_g, chr_b[i_b], end_b[i] + window + 1);
				}
				block_id++;
			}
		}
		i_b = i_b_e + 1;
		// std::cerr << block_id << " blocks ...\n";
	}
}


