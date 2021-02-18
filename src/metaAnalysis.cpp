/*  
    Copyright (C) 2020 
    Authors: Corbin Quick <qcorbin@hsph.harvard.edu>
	         Li Guan <guanli@umich.edu>

    This file is a part of APEX.

    APEX is distributed "AS IS" in the hope that it will be 
    useful, but WITHOUT ANY WARRANTY; without even the implied 
    warranty of MERCHANTABILITY, NON-INFRINGEMENT, or FITNESS 
    FOR A PARTICULAR PURPOSE.

    The above copyright notice and disclaimer of warranty must 
    be included in all copies or substantial portions of APEX.
*/

#include "metaAnalysis.hpp"


void vcov_meta_data::open(const std::vector<std::string>& study_prefixes, const std::string& region )
{
	int i = 0;
	for( const std::string& prefix : study_prefixes ){
		i++;
		std::cerr << "Processing Study #" << i << " ...";
		vc.push_back(
			vcov_data(prefix, region)
		);
		clear_line_cerr();
		std::cerr << "Processed Study #" << i << " (" << prefix << ") with " << vc[i-1].chr.size() << " variants.\n";
		
	}
	std::cerr << "\n";
	build_index();
}

void vcov_meta_data::build_index()
{
	
	int n_studies = vc.size();
	
	std::vector<int> ii(n_studies, 0);
	
	std::vector<int> n_var(n_studies);
	
	si.clear(); 
	sflipped.clear();
	
	for( int i = 0; i < n_studies; i++ ){
		si.push_back(std::vector<int>(0));
		sflipped.push_back(std::vector<bool>(0));
		n_var[i] = vc[i].chr.size();
		
		//cerr << "Study #" << i << ": " << n_var[i] << " variants.\n";
	}
	std::string iter_cerr_suffix = " shared variants ...";
	
	std::cerr << "Found ";
	print_iter_cerr(1, 0, iter_cerr_suffix);
	int last = 0;
	
	while( all_lt( ii, n_var ) ){

		bool skip = false;
		
		int pos_0 = 0;
		int ichr_0 = i_chrom( vc[0].chr[ii[0]] ); 
		
		// find the maximum current chromosome
		for( int s = 0; s < n_studies; s++ ){
			if( i_chrom( vc[s].chr[ii[s]] ) > ichr_0 ){
				ichr_0 = i_chrom( vc[s].chr[ii[s]] );
			}
		}
		
		// seq all studies to the max chromosome
		for( int s = 0; s < n_studies; s++ ){
			while( i_chrom( vc[s].chr[ii[s]] ) < ichr_0 ){
				ii[s]++;
				if( ii[s] >= n_var[s] ) break;
			}
			if( i_chrom( vc[s].chr[ii[s]] ) != ichr_0 ){
				skip = true;
			}
		}
		
		// break if any study reached the end
		if( !all_lt( ii, n_var ) ) break;
		// skip to next candidate if the chromosome isn't shared across studies
		if( skip ) continue;
		
		// find the maximum current position
		for( int s = 0; s < n_studies; s++ ){
			if( vc[s].pos[ii[s]] > pos_0 ){
				pos_0 = vc[s].pos[ii[s]];
			}
		}
		
		
		bool multi_allelic = false;
		
		// seq all studies to the max position
		for( int s = 0; s < n_studies; s++ ){
			while( vc[s].pos[ii[s]] < pos_0 ){
				ii[s]++;
				if( ii[s] >= n_var[s] ) break;
				if( i_chrom( vc[s].chr[ii[s]] ) != ichr_0 ){
					skip = true;
					// break because we're on a new chromosome
					break;
				}
			}
			if( vc[s].pos[ii[s]] != pos_0 ){
				// if the position is absent in a single study, we skip it 
				skip = true;
				// we continue the outer loop to move all study indexes forward
			}else if( global_opts::biallelic_only && ii[s] + 1 < vc[s].pos.size() && !multi_allelic ){
				// check if the next variant has the same position. 
				// if so, it is a multi-allelic site 
				multi_allelic = (
					(vc[s].pos[ii[s] + 1] == pos_0) && 
					(i_chrom( vc[s].chr[ii[s]+ 1 ]) == ichr_0)
				);
			}
		}
		
		// break if any study reached the end
		if( !all_lt( ii, n_var ) ) break;
		// skip to next candidate if the position isn't shared across studies
		if( skip ) continue;
		
		// OK, we can now be sure that all studies have the same chr,pos.
		
		if( global_opts::biallelic_only && multi_allelic ){
			// if we're only keeping biallelic variants,
			// and the current site has >2 alleles, move 
			// all studies forward and then continue. 
	
			// seq all studies to the next position.
			for( int s = 0; s < n_studies; s++ ){
				while( vc[s].pos[ii[s]] <= pos_0 && i_chrom(vc[s].chr[ii[s]])==ichr_0 ){
					ii[s]++;
					if( ii[s] >= n_var[s] ) break;
				}
			}
			// back to start.
			continue; 
		}
		
		std::string chr_0 = vc[0].chr[ii[0]];
		std::string alt_0 = vc[0].alt[ii[0]];
		std::string ref_0 = vc[0].ref[ii[0]];
		
		for( int s = 0; s < n_studies; s++ ){
			if( !( vc[s].ref[ii[s]] == ref_0 && vc[s].alt[ii[s]] == alt_0  ) ){
				ii[s]++;
				skip = true;
			}
		}
		// skip to next candidate if the alleles don't match across studies
		if( skip ) continue;
		
		// OK, now we're sure that all variants match
		chr.push_back(chr_0);
		pos.push_back(pos_0);
		ref.push_back(ref_0);
		alt.push_back(alt_0);
		
		double var_0 = 0.0, mac_0 = 0.0;
		std::vector<double> var_perStudy_0(n_studies);
		
		for( int s = 0; s < n_studies; s++ ){
		
			// location of the shared variant in study s's variant list
			si[s].push_back(ii[s]);
			
			int ns = si[s].size();
			if( ns > 1 ){
				if( si[s][ns-1] < si[s][ns-2] ){
					std::cout << "FATAL\n";
					std::cout << s << ": " << ns << ", " << si[s][ns-1] << ", " << si[s][ns-2] << "\n";
					abort();
				}
			}
			
			var_perStudy_0[s] = vc[s].var[ii[s]];
			
			// is it flipped relative to study #1? (Currently ignored; we match ref & alt)
			sflipped[s].push_back(false);
			
			mac_0 += vc[s].mac[ii[s]];
			var_0 += vc[s].var[ii[s]];
			
			// increment forward by 1
			ii[s]++;
		}
		var.push_back(var_0);
		var_perStudy.push_back(var_perStudy_0);
		mac.push_back(mac_0);
		
		thinned_iter_cerr(last, chr.size(), iter_cerr_suffix);
	}
	clear_line_cerr();
	std::cerr << "Meta-analysis: " << chr.size() << " total shared variants.\n";
}

Eigen::MatrixXd vcov_meta_data::getV_perStudy(const int& s, const std::vector<int>& v_i)
{
	int n = v_i.size();
	std::vector<int> v_s(n);
	
	for( int j = 0; j < n; ++j){
		v_s[j] = si[s][v_i[j]];
	}
	return vc[s].getV(v_s);
}


Eigen::MatrixXd vcov_meta_data::getGtG_perStudy(const int& s, const std::vector<int>& v_i, const std::vector<int>& v_j)
{
	int n = v_i.size();
	int m = v_j.size();

	std::vector<int> v_i_s(n);
	std::vector<int> v_j_s(m);
	
	for( int k = 0; k < n; ++k){
		v_i_s[k] = si[s][v_i[k]];
	}
	for( int k = 0; k < m; ++k){
		v_j_s[k] = si[s][v_j[k]];
	}
	return vc[s].getGtG(v_i_s, v_j_s);
}

Eigen::MatrixXd vcov_meta_data::getV_perStudy(const int& s, const std::vector<int>& v_i, const std::vector<int>& v_j)
{
	int n = v_i.size();
	int m = v_j.size();

	std::vector<int> v_i_s(n);
	std::vector<int> v_j_s(m);
	
	for( int k = 0; k < n; ++k){
		v_i_s[k] = si[s][v_i[k]];
	}
	for( int k = 0; k < m; ++k){
		v_j_s[k] = si[s][v_j[k]];
	}
	return vc[s].getV(v_i_s, v_j_s);
}

Eigen::MatrixXd vcov_meta_data::getGtG_perStudy(const int& s, const std::vector<int>& v_i)
{
	int n = v_i.size();
	std::vector<int> v_s(n);
	
	for( int j = 0; j < n; j++){
		v_s[j] = si[s][v_i[j]];
	}
	return vc[s].getGtG(v_s);
}

Eigen::MatrixXd vcov_meta_data::getGtG(const std::vector<int>& v_i, const std::vector<double>& w, const std::vector<int>& sl )
{
	int n = v_i.size();
	Eigen::MatrixXd out = Eigen::MatrixXd::Zero(n, n);

	if( sl.size() == 0){
		for( int s = 0; s < vc.size(); s++)
		{
			double w_s = (w.size() == 0) ? 1 : w[s];
			out += w_s*getGtG_perStudy(s, v_i);
		}	
	}else{
		for( const int& s : sl)
		{
			double w_s = (w.size() == 0) ? 1 : w[s];
			out += w_s*getGtG_perStudy(s, v_i);
		}	
	}
	
	return out;
}

Eigen::MatrixXd vcov_meta_data::getV(const std::vector<int>& v_i, const std::vector<double>& w, const std::vector<int>& sl )
{
	int n = v_i.size();
	Eigen::MatrixXd out = Eigen::MatrixXd::Zero(n, n);

	if( sl.size() == 0){
		for( int s = 0; s < vc.size(); s++)
		{
			double w_s = (w.size() == 0) ? 1 : w[s];
			out += w_s * getV_perStudy(s, v_i);
		}
	}else{
		for( const int& s : sl )
		{
			double w_s = (w.size() == 0) ? 1 : w[s];
			out += w_s * getV_perStudy(s, v_i);
		}
	}
	
	return out;
}

Eigen::MatrixXd vcov_meta_data::getGtG(const std::vector<int>& v_i, const std::vector<int>& v_j, const std::vector<double>& w, const std::vector<int>& sl )
{
	int n = v_i.size();
	int m = v_j.size();
	Eigen::MatrixXd out = Eigen::MatrixXd::Zero(n, m);

	if( sl.size() == 0){
		for( int s = 0; s < vc.size(); s++)
		{
			double w_s = (w.size() == 0) ? 1 : w[s];
			out += w_s * getGtG_perStudy(s, v_i, v_j);
		}
	}else{
		for( const int& s : sl )
		{
			double w_s = (w.size() == 0) ? 1 : w[s];
			out += w_s * getGtG_perStudy(s, v_i, v_j);
		}
	}
	
	return out;
}

Eigen::MatrixXd vcov_meta_data::getV(const std::vector<int>& v_i, const std::vector<int>& v_j, const std::vector<double>& w, const std::vector<int>& sl )
{
	int n = v_i.size();
	int m = v_j.size();
	Eigen::MatrixXd out = Eigen::MatrixXd::Zero(n, m);
	
	if( sl.size() == 0){
		for( int s = 0; s < vc.size(); s++)
		{
			double w_s = (w.size() == 0) ? 1 : w[s];
			out += w_s * getV_perStudy(s, v_i, v_j);
		}
	}else{
		for( const int& s : sl )
		{
			double w_s = (w.size() == 0) ? 1 : w[s];
			out += w_s * getV_perStudy(s, v_i, v_j);
		}
	}
	
	return out;
}

cis_sumstat_data::cis_sumstat_data(const std::string& pf, const std::string& reg)
{
	open(pf, reg);
}

void format_gene_ids(std::vector<std::string>& ids)
{
	if( global_opts::trim_gene_ids )
	{
		for(std::string& id : ids){
			id = id.substr(0, id.find("."));
		}
	}	
}

void cis_sumstat_data::open(const std::string& pf, const std::string& reg)
{
	file_prefix = pf;
	region = reg;
	
	std::string sumstat_file = file_prefix + ".cis_sumstats.txt.gz";
	std::string gene_file = file_prefix + ".cis_gene_table.txt.gz";
	
	kstring_t str = {0,0,0};

	int ncols = 0;
	n_genes = 0;
	
	data_parser dp;
	
	dp.add_field(chr,0);
	dp.add_field(start,1);
	dp.add_field(end,2);
	dp.add_field(gene_id,3);
	dp.add_field(egene_pval,4);
	dp.add_field(NS,5);
	dp.add_field(NC,6);
	dp.add_field(SD,7);
	dp.add_field(N_CIS,8);
	
	// dp.parse_file(gene_file, region);
	
	indexed_hts_file htsf(gene_file, region);
	int nrows = 0;
	dp.parse_file(htsf, nrows);
	htsf.close();
	
	format_gene_ids(gene_id);
	
	dp.clear();
	
	// dp.add_field(start,1);
	// dp.add_field(end,2);
	dp.add_field(S_CIS,3);	
	
	indexed_hts_file ss_hfile(sumstat_file, region);
	
	int n_matched = 0;
	
	while( ss_hfile.next_line(str) >= 0 && n_genes < N_CIS.size() )
	{
		if( !str.l ) break;
		if ( str.s[0] == '#' ) continue;
		
		int n_fields;
	
		int *offsets = ksplit(&str, 0, &n_fields);
	
		if( N_CIS[n_genes] != n_fields - 4 ){
			if( n_matched > 0 ){
				std::cerr << "MISMATCHED: " << N_CIS[n_genes] << ", " << n_fields - 4 << "\n";
				abort();
			}
			continue;
		}else{
			n_matched++;
		}
		
		dp.parse_fields(str, offsets, n_fields);
		
		// std::cout << S_CIS[n_genes] << "\n";
		
		/*if( global_opts::filter_genes ){
			if(!has_element(global_opts::target_genes, gene_id[n_genes])){
				n_genes++;
				continue;
			}
		}*/
		
		Eigen::VectorXd score_v(N_CIS[n_genes]);
		
		for( int i = 0; i < N_CIS[n_genes]; ++i )
		{
			score_v(i) = atof(str.s + offsets[i + 4]);
		}
		
		score.push_back(score_v);
		n_genes++;
	}
	
	ks_free(&str);
	ss_hfile.close();
	
	std::cerr << "Processed cis summary statistics for " <<  chr.size() << " genes.\n";
	
	// ln.set(vcov.pos);
}

void cis_meta_data::merge_intersection(const std::vector<std::vector<int>>& si, const std::vector<std::vector<bool>>& sflipped)
{
	// Merge sumstat data, keeping the intersection of genes across studies.
	
	int n_studies = ss.size();
	
	// gene index
	std::vector<int> jj(n_studies, 0);
	
	// number of genes per study
	std::vector<int> n_genes(n_studies);
	
	// start and end indexes
	int i_s = 0, i_e = 0;
	
	for( int i = 0; i < n_studies; i++ ){
		n_genes[i] = ss[i].gene_id.size();
	}
	std::string iter_cerr_suffix = " shared genes ...";
	
	std::cerr << "Found ";
	print_iter_cerr(1, 0, iter_cerr_suffix);
	int last = 0;
	
	while( all_lt( jj, n_genes ) ){

		bool skip = false;
		
		int start_0 = 0;
		int end_0 = 0;
		int ichr_0 = i_chrom( ss[0].chr[jj[0]] ); 
		
		// find the maximum current chromosome
		for( int s = 0; s < n_studies; s++ ){
			if( i_chrom( ss[s].chr[jj[s]] ) > ichr_0 ){
				ichr_0 = i_chrom( ss[s].chr[jj[s]] );
			}
		}
		
		// seq all studies to the max chromosome
		for( int s = 0; s < n_studies; s++ ){
			while( i_chrom( ss[s].chr[jj[s]] ) < ichr_0 ){
				jj[s]++;
				if( jj[s] >= n_genes[s] ) break;
			}
			if( i_chrom( ss[s].chr[jj[s]] ) != ichr_0 ){
				skip = true;
			}
		}
		
		// break if any study reached the end
		if( !all_lt( jj, n_genes ) ) break;
		// skip to next candidate if the chromosome isn't shared across studies
		if( skip ) continue;
		
		// find the maximum current position
		for( int s = 0; s < n_studies; s++ ){
			if( ss[s].start[jj[s]] > start_0 ){
				start_0 = ss[s].start[jj[s]];
				end_0 = ss[s].end[jj[s]];
			}
		}
		
		// std::cout << start_0 << "\n";
		
		// seq all studies to the max position
		for( int s = 0; s < n_studies; s++ ){
			while( ss[s].start[jj[s]] < start_0 ){
				jj[s]++;
				if( jj[s] >= n_genes[s] ) break;
				if( i_chrom( ss[s].chr[jj[s]] ) != ichr_0 ){
					skip = true;
					// break because we're on a new chromosome
					break;
				}
			}
			if( ss[s].start[jj[s]] != start_0 ){
				// if the position is absent in a single study, we skip it 
				skip = true;
				// we continue the outer loop to move all study indexes forward
			}
		}
		// break if any study reached the end
		if( !all_lt( jj, n_genes ) ) break;
		// skip to next candidate if the position isn't shared across studies
		if( skip ) continue;
		
		// OK, we can now be sure that all studies have the same chr,start.
		
		std::string chr_0 = ss[0].chr[jj[0]];
		std::string gene_id_0 = ss[0].gene_id[jj[0]];
		
		for( int s = 0; s < n_studies; s++ ){
			if( !( ss[s].gene_id[jj[s]] == gene_id_0 ) ){
				jj[s]++;
				skip = true;
			}
		}
		
		/*if( !skip && global_opts::filter_genes ){
			if(!has_element(global_opts::target_genes, gene_id_0)){
				for( int s = 0; s < n_studies; s++ ){
					jj[s]++;
				}
				skip = true;
			}
		}*/
		
		// skip to next candidate if the gene name doesn't match across studies
		
		if( skip ) continue;
		
		double N_0 = 0;
		double DF_0 = 0;
		double SD_0 = 0; 
		double DENOM = 0;
		
		std::vector<double> SD_perStudy_0(n_studies);
		std::vector<double> DF_perStudy_0(n_studies);
		std::vector<double> SSR_perStudy_0(n_studies);
		
		std::vector<double> ivw_0(n_studies);
		
		for( int s = 0; s < n_studies; s++ ){
			
			SD_perStudy_0[s] = ss[s].SD[jj[s]];
			DF_perStudy_0[s] = ss[s].NS[jj[s]] - ss[s].NC[jj[s]];
			
			SSR_perStudy_0[s] = (ss[s].NS[jj[s]] - 1)*SD_perStudy_0[s]*SD_perStudy_0[s];
			
			ivw_0[s] = DF_perStudy_0[s]/SSR_perStudy_0[s];
			
			N_0 += ss[s].NS[jj[s]];
			DF_0 += DF_perStudy_0[s];
			SD_0 += (ss[s].NS[jj[s]] - 1)*SD_perStudy_0[s]*SD_perStudy_0[s] * ivw_0[s];
			DENOM += ivw_0[s];
			
		}
		SD_0 = std::sqrt(SD_0 / DENOM);

		std::vector<int> c_s(n_studies);
		std::vector<int> c_e(n_studies);
		
		
		// now we want to find the range of overlapping indexes
		for( int s = 0; s < n_studies; s++ ){
			c_s[s] = ss[s].S_CIS[jj[s]];
			c_e[s] = ss[s].S_CIS[jj[s]] + ss[s].N_CIS[jj[s]] - 1;
			while( si[s][i_s] < c_s[s] && i_s < si[s].size() ){
				i_s++;
			}
			while( si[s][i_e] < c_e[s] && i_e < si[s].size() ){
				i_e++;
			}			
		}
		// shift the end index back, so that we have the intersection of variants
		for( int s = 0; s < n_studies; s++ ){
			while( si[s][i_e] > c_e[s] || i_e >= si[s].size() ){
				i_e--;
			}			
		}

		Eigen::VectorXd score_0 = Eigen::VectorXd::Zero(i_e - i_s + 1);
		std::vector<Eigen::VectorXd> score_perStudy_0(n_studies, Eigen::VectorXd::Zero(i_e - i_s + 1));
		//Eigen::VectorXd var_0 = Eigen::VectorXd::Zero(i_e - i_s + 1);
		
		for( int s = 0; s < n_studies; s++ ){
			const Eigen::VectorXd& sc_s = ss[s].score[jj[s]];
			for( int i = 0, ii = i_s; ii <= i_e; ++ii, ++i ){
				int idx = si[s][ii] - ss[s].S_CIS[jj[s]];
				
				if( idx < 0 || idx >= sc_s.size() ){
					std::cerr << "\nFatal: index out of bounds in cis_meta_data::merge\n";
					abort();
				}
				
				double adj_s = ivw_0[s]*SD_perStudy_0[s];
				
				if( sflipped[s][ii] ){
					score_0(i) -= sc_s( idx )*adj_s;
					score_perStudy_0[s](i) -= sc_s( idx )*SD_perStudy_0[s];
				}else{
					score_0(i) += sc_s( idx )*adj_s;
					score_perStudy_0[s](i) += sc_s( idx )*SD_perStudy_0[s];
				}
			}
			// now increment to the next gene
			jj[s]++;
		}
		// move forward if there are no shared variants
		// if( score_0.size() == 0 ) continue;
		
		// std::cout << start_0 << "\t" << gene_id_0 << "\n";

		// OK, now we're sure that the gene is shared across studies
		
		/*
		if( global_opts::filter_genes ){
			if( !has_element(global_opts::target_genes, gene_id_0) ){
				continue;
			}
		}
		*/
		
		chr.push_back(chr_0);
		start.push_back(start_0);
		end.push_back(end_0);
		gene_id.push_back(gene_id_0);
		
		S_CIS.push_back(i_s);
		N_CIS.push_back(i_e - i_s + 1);
		
		N.push_back(N_0);
		DF.push_back(DF_0);
		SD.push_back(SD_0);
		score.push_back(score_0);
		score_perStudy.push_back(score_perStudy_0);
		
		SSR_perStudy.push_back(SSR_perStudy_0);
		SD_perStudy.push_back(SD_perStudy_0);
		DF_perStudy.push_back(DF_perStudy_0);
		ivw.push_back(ivw_0);
		
		std::vector<int> all_studies;
		for(int i = 0; i < n_studies; i++){
			all_studies.push_back(i);
		}
		study_list.push_back(all_studies);
		
		// var_score.push_back(var_0);
		
		thinned_iter_cerr(last, chr.size(), iter_cerr_suffix, 1);
	}
	
	clear_line_cerr();
	std::cerr << "Meta-analysis: " << chr.size() << " total shared genes.\n";
}


void cis_meta_data::merge(const std::vector<std::vector<int>>& si, const std::vector<std::vector<bool>>& sflipped)
{
	// Merge sumstat data, keeping the union of genes across studies.
	
	int n_studies = ss.size();
	
	// gene index
	std::vector<int> jj(n_studies, 0);
	
	// number of genes per study
	std::vector<int> n_genes(n_studies);
	
	// start and end indexes
	int i_s = 0, i_e = 0;
	
	for( int i = 0; i < n_studies; i++ ){
		n_genes[i] = ss[i].gene_id.size();
	}
	
	std::string iter_cerr_suffix = " total genes ...";
	
	std::cerr << "Found ";
	print_iter_cerr(1, 0, iter_cerr_suffix);
	int last = 0;
	
	std::vector<std::string> all_gene_id;
	
	std::vector<int> c_studies = which_lt(jj, n_genes);
	
	while( c_studies.size() > 0 ){
		
		int c_0 = c_studies[0];
		
		std::vector<int> move_back(c_studies.size(), 0);
		
		int start_0 = ss[c_0].start[jj[c_0]];
		int end_0 = ss[c_0].end[jj[c_0]];
		int ichr_0 = i_chrom(ss[c_0].chr[jj[c_0]]);
		
		int s_ref = c_0;
		
		// find the minimum current chromosome and position.
		for( const int& s : c_studies ){
			if( i_chrom( ss[s].chr[jj[s]] ) < ichr_0 ){
				ichr_0 = i_chrom( ss[s].chr[jj[s]] );
				start_0 = ss[s].start[jj[s]];
				end_0 = ss[s].end[jj[s]];
				s_ref = s;
			}else if( i_chrom( ss[s].chr[jj[s]] ) == ichr_0 && ss[s].start[jj[s]] <= start_0 ){
				start_0 = ss[s].start[jj[s]];
				end_0 = ss[s].end[jj[s]];
				s_ref = s;
			}
		}
		
		std::string chr_0 = ss[s_ref].chr[jj[s_ref]];
		std::string gene_id_0 = ss[s_ref].gene_id[jj[s_ref]];
		
		std::vector<int> studies_with_gene;
		
		for( const int& s : c_studies ){
			if( ss[s].gene_id[jj[s]] == gene_id_0 ){
				studies_with_gene.push_back(s);
			}else{
				while( ss[s].gene_id[jj[s]] != gene_id_0 && jj[s] > 0 && ss[s].start[jj[s]] == start_0 ){
					if( ss[s].start[jj[s] - 1] == start_0 ){
						jj[s]--;
						move_back[s]--;
						if( ss[s].gene_id[jj[s]] == gene_id_0 ){
							studies_with_gene.push_back(s);
						}
					}else{
						break;
					}
				}
				while( ss[s].gene_id[jj[s]] != gene_id_0 && jj[s] < n_genes[s] && ss[s].start[jj[s]] == start_0 ){
					if( ss[s].start[jj[s] + 1] == start_0 ){
						jj[s]++;
						move_back[s]++;
						if( ss[s].gene_id[jj[s]] == gene_id_0 ){
							studies_with_gene.push_back(s);
						}
					}else{
						break;
					}
				}
			}
		}
		
		double N_0 = 0;
		double DF_0 = 0;
		double SD_0 = 0; 
		double DENOM = 0;
		
		std::vector<double> SD_perStudy_0(n_studies);
		std::vector<double> DF_perStudy_0(n_studies);
		std::vector<double> SSR_perStudy_0(n_studies);
		
		std::vector<double> ivw_0(n_studies);
		
		for( const int& s : studies_with_gene ){
			
			SD_perStudy_0[s] = ss[s].SD[jj[s]];
			DF_perStudy_0[s] = ss[s].NS[jj[s]] - ss[s].NC[jj[s]];
			
			SSR_perStudy_0[s] = (ss[s].NS[jj[s]] - 1)*SD_perStudy_0[s]*SD_perStudy_0[s];
			
			ivw_0[s] = DF_perStudy_0[s]/SSR_perStudy_0[s];
			
			N_0 += ss[s].NS[jj[s]];
			DF_0 += DF_perStudy_0[s];
			SD_0 += (ss[s].NS[jj[s]] - 1)*SD_perStudy_0[s]*SD_perStudy_0[s] * ivw_0[s];
			DENOM += ivw_0[s];
			
		}
		SD_0 = std::sqrt(SD_0 / DENOM);

		std::vector<int> c_s(n_studies);
		std::vector<int> c_e(n_studies);
		
		// std::cout << "\nStudies with gene " << gene_id_0 <<  " : " << studies_with_gene.size() << "\n";
		
		// move back indices if needed
		for( const int& s : studies_with_gene ){
			// std::cout << s << " : jj[s] = " << jj[s] << "\n";
			c_s[s] = ss[s].S_CIS[jj[s]];
			c_e[s] = ss[s].S_CIS[jj[s]] + ss[s].N_CIS[jj[s]] - 1;
			while( si[s][i_s] >= c_s[s] && i_s > 0 ){
				i_s--;
			}
			while( si[s][i_e] < c_e[s] && i_e < si[s].size() ){
				i_e++;
			}			
		}
		
		// std::cout << i_s << " " << i_e << "\n";
		
		// now we want to find the range of overlapping indexes
		for( const int& s : studies_with_gene ){
			while( i_s >= si[s].size() && i_s > 0){
				i_s--;
			}
			while( si[s][i_s] < c_s[s] && i_s < si[s].size() ){
				i_s++;
			}
			while( i_e >= si[s].size() && i_e > 0){
				i_e--;
			}
			while( si[s][i_e] >= c_e[s] && i_e > 0 ){
				i_e--;
			}
		}
		
		// std::cout << i_s << " " << i_e << "\n";
		
		
		bool has_zero_snps = false;
		if( i_e - i_s < 0 || i_e + i_s == 0 ){
			has_zero_snps = true;
			std::cerr << "\nWarning: Gene " << gene_id_0 << " has 0 variants after merging.\n";
		}
		
		Eigen::VectorXd score_0 = Eigen::VectorXd::Zero(i_e - i_s + 1);
		std::vector<Eigen::VectorXd> score_perStudy_0(n_studies, Eigen::VectorXd::Zero(i_e - i_s + 1));
		//Eigen::VectorXd var_0 = Eigen::VectorXd::Zero(i_e - i_s + 1);
		
		
		
		for( const int& s : studies_with_gene ){
			
			
			// std::cout << i_s << "\t" << i_e << "\n";
			// std::cout << jj[s] << "\n";
			
			const Eigen::VectorXd& sc_s = ss[s].score[jj[s]];
		
			if( !has_zero_snps ){
				for( int i = 0, ii = i_s; ii <= i_e; ii++, i++ ){
					
					int idx = si[s][ii] - ss[s].S_CIS[jj[s]];
					
					if( idx < 0 || idx >= sc_s.size() ){
						std::cerr << "\nFatal: index out of bounds in cis_meta_data::merge\n";
						abort();
					}
					
					double adj_s = ivw_0[s]*SD_perStudy_0[s];
					
					if( sflipped[s][ii] ){
						score_0(i) -= sc_s( idx )*adj_s;
						score_perStudy_0[s](i) -= sc_s( idx )*SD_perStudy_0[s];
					}else{
						score_0(i) += sc_s( idx )*adj_s;
						score_perStudy_0[s](i) += sc_s( idx )*SD_perStudy_0[s];
					}
				}
			}
			// now increment to the next gene
			if( move_back[s] == 0 ){
				jj[s]++;
			}
		}
		for( const int& s : c_studies ){
			jj[s] -= move_back[s];
			if( jj[s] < 0 ){
				jj[s] = 0;
			}
		}
		
		all_gene_id.push_back(gene_id_0);
		
		for( const int& s : c_studies ){
			while( jj[s] < n_genes[s] ){
				if( has_element(all_gene_id, ss[s].gene_id[jj[s]]) ){
					jj[s]++;
				}else{
					break;
				}
			}
		}
		
		if( has_zero_snps ){
			continue;
		}
		
		/*
		if( global_opts::filter_genes ){
			if( !has_element(global_opts::target_genes, gene_id_0) ){
				continue;
			}
		}
		*/
		
		chr.push_back(chr_0);
		start.push_back(start_0);
		end.push_back(end_0);
		gene_id.push_back(gene_id_0);
		
		S_CIS.push_back(i_s);
		N_CIS.push_back(i_e - i_s + 1);
		
		N.push_back(N_0);
		DF.push_back(DF_0);
		SD.push_back(SD_0);
		score.push_back(score_0);
		score_perStudy.push_back(score_perStudy_0);
		
		SSR_perStudy.push_back(SSR_perStudy_0);
		SD_perStudy.push_back(SD_perStudy_0);
		DF_perStudy.push_back(DF_perStudy_0);
		ivw.push_back(ivw_0);
		
		study_list.push_back(studies_with_gene);
		
		// var_score.push_back(var_0);
		
		c_studies = which_lt(jj, n_genes);
		
		thinned_iter_cerr(last, chr.size(), iter_cerr_suffix, 1);
	}
	
	clear_line_cerr();
	std::cerr << "Meta-analysis: " << chr.size() << " total genes across studies.\n";
}


void cis_meta_data::meta_analyze()
{
	std::string out_name = global_opts::out_prefix + ".cis_meta.svar.tsv";
	std::ofstream os(out_name.c_str(), std::ofstream::out);

	std::vector<std::string> col_names{"#chr", "pos", "ref", "alt", "gene", "studies", "beta", "se", "pval"};
	
	print_header(col_names, os);
	
	for(int i = 0; i < chr.size(); ++i){
		std::string in_studies = std::to_string(study_list[i][0] + 1);
		for( int s = 1; s < study_list[i].size(); s++){
			in_studies += ("," + std::to_string(study_list[i][s] + 1));
		}
		
		for(int j = S_CIS[i], jj = 0; jj < N_CIS[i]; ++j, ++jj ){
			
			double beta, se;
			
			if( global_opts::IVW_H1_SIGMA ){
				beta = 0;
				double denom = 0;
				for( const int& s : study_list[i] ){
					
					double sc = score_perStudy[i][s](jj);
					double dv = diagV_perStudy(i,jj,s);
					
					double sigma_sq = (SSR_perStudy[i][s] - sc*sc/dv)/(DF_perStudy[i][s]-1);
					
					denom += dv/sigma_sq;
					beta += sc/sigma_sq;
				}
				beta = beta/denom;
				se = std::sqrt(1/denom);
				
			}else{
				double dv = diagV(i, jj);
				if(dv > 0 ){
					beta = score[i](jj)/dv;
					se = std::sqrt( (DF[i] - 1)/dv - beta*beta)/std::sqrt(DF[i] - 1);
				}
			}

			if( se > 0  ){

				double z = beta/se;
				double pval = pf( z*z, 1.0, DF[i] - 1, true );
				
				os << 
					//score_perStudy[i][0](jj) << ", " << score_perStudy[i][1](jj) << "\t" <<
					//SD[i] << " " << dv << " " << (DF[i] - 1) << "\t" <<
					chr[i] << "\t" <<
					pos(j) << "\t" <<
					ref(j) << "\t" <<
					alt(j) << "\t" <<
					gene_id[i] << "\t" <<
					in_studies << "\t" << 
	//				SD[i]*beta << "\t" << 
	//				SD[i]*se << "\t" << 
					beta << "\t" << 
					se << "\t" << 
					pval << "\n";
			}
		}
	}
	
	os.close();
}

std::vector<int> seq_from(const int& s, const int& n)
{
	std::vector<int> out(n);
	for( int i = 0; i < n; i++){
		out[i] = s + i;
	}
	return out;
}

const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ",", "\n");

void cis_meta_data::get_vcov_gene(const int& gene_index, const bool& centered)
{
	
	if( N_CIS[gene_index] <= 0 ){
		std::cerr << "\nERROR: No variants present for " << gene_id[gene_index] << ".\n";
		std::cerr << "Note: score.size() == " << score[gene_index].size() << "\n\n";
		
		return;
	}
	
	//cout << "0\n";
	
	std::string& gene = gene_id[gene_index];
	Eigen::VectorXd U = score[gene_index];
	
	double df0 = DF[gene_index];
	double n = N[gene_index];
	
	double stdev = SD[gene_index];
	
	int n_var = N_CIS[gene_index];
	int s_var = S_CIS[gene_index];
	
	Eigen::VectorXd dV( n_var );
	
	//cout << "1\n";
	
	if( U.size() != n_var ){
		std::cerr << "U.size() != n_var\n";
		std::cerr << U.size() << " != " << n_var << "\n";
		std::cerr << U(0) << "\t" << U(U.size()-1) << "\n";
		exit(1);
	}
	
	for( int i = 0; i < n_var; i++ ){
		//dV(i) = vc.var[i + s_var];
		dV(i) = diagV(gene_index, i);
	}
	
	vcov_getter vget(vc, ivw[gene_index], s_var, n_var, study_list[gene_index]);
	
	auto snp = [&](const int& i ){ int j = s_var + i; return vc.chr[j] + "_" + std::to_string(vc.pos[j]) + "_" + vc.ref[j] + "_" + vc.alt[j];};
	
	for(int i = 0; i < n_var; i++){
		if( i + 1 < n_var){
			std::cout << snp(i) << ",";
		}else{
			std::cout << snp(i) << "\n";
		}
	}
	
	if( centered ){
		std::cout << vget.Var(seq_from(0,n_var)).format(CSVFormat) << "\n";
	}else{
		std::cout << vget.Var_uncentered(seq_from(0,n_var)).format(CSVFormat) << "\n";
	}
	
	return;
}


void cis_meta_data::conditional_analysis(const int& gene_index, std::ostream& os, std::ostream& os_b, std::ostream& os_log )
{
	
	if( N_CIS[gene_index] <= 0 ){
		std::cerr << "\nWARNING: No variants present for " << gene_id[gene_index] << ".\n";
		std::cerr << "    Note: score.size() == " << score[gene_index].size() << "\n\n";
		
		os_log << "\nWARNING: No variants present for " << gene_id[gene_index] << ".\n";
		os_log << "    Note: score.size() == " << score[gene_index].size() << "\n\n";

		return;
	}
	
	//cout << "0\n";
	
	std::string& gene = gene_id[gene_index];
	Eigen::VectorXd U = score[gene_index];
	
	double df0 = DF[gene_index];
	double n = N[gene_index];
	
	double stdev = SD[gene_index];
	
	int n_var = N_CIS[gene_index];
	int s_var = S_CIS[gene_index];
	
	// distance to TSS
	double tss_pos = 0.5*(start[gene_index] + end[gene_index]);
	std::vector<double> dtss(n_var);
	for(int i = 0; i < n_var; i++){
		dtss[i] = ( tss_pos - (double) vc.pos[s_var + i]);
	}
	
	Eigen::VectorXd dV( n_var );
	
	if( U.size() != n_var ){
		std::cerr << "U.size() != n_var\n";
		std::cerr << U.size() << " != " << n_var << "\n";
		std::cerr << U(0) << "\t" << U(U.size()-1) << "\n";
		exit(1);
	}
	
	for( int i = 0; i < n_var; i++ ){
		//dV(i) = vc.var[i + s_var];
		dV(i) = diagV(gene_index, i);
	}
	
	vcov_getter vget(vc, ivw[gene_index], s_var, n_var, study_list[gene_index]);
	
	stdev = 1;

	forward_lm out(U, dV, n, df0, stdev, vget, global_opts::LM_ALPHA, dtss);
	
	//cout << "4\n";
	std::string in_studies = std::to_string(study_list[gene_index][0] + 1);
	for( int s = 1; s < study_list[gene_index].size(); s++){
		in_studies += ("," + std::to_string(study_list[gene_index][s] + 1));
	}

	if( out.error_message != "" ){
		os_log << "ERROR: GENE = " + gene + " :\n" + out.error_message + "\n"; 
	}

	auto snp = [&](const int& i ){ int j = s_var + i; return vc.chr[j] + "_" + std::to_string(vc.pos[j]) + "_" + vc.ref[j] + "_" + vc.alt[j];};
	
	if( out.beta.size() > 0 ){
		
		for(int i = 0; i < out.beta.size(); ++i)
		{
			// os.precision(4);
			os << gene << "\t" << in_studies;
			os << "\t" << i+1 << ":" << out.beta.size() << "\t" << snp(out.keep[i]) << "\t" << 
			out.beta[i] << "\t" << out.se[i] << "\t" << out.pval_joint[i]  << "\t" << out.pval_adj[i] << "\t";
			// os.precision(2);
			os << out.pval_0[i] << "\t" << out.pval_seq[i] << "\n";
		}
		// os.precision(5);
		
		if( global_opts::RSQ_BUDDY < 1.00 ){
			for(int i = 0; i < out.beta.size(); ++i){
				for(int j = 0; j < out.buddy_list[i].size(); j++){
					os_b << 
						snp(out.keep[i]) << "\t" << 
						snp(out.buddy_list[i][j]) << "\t" << 
						out.corr_list[i][j] << "\t" << 
						out.corr_list[i][j] * out.corr_list[i][j] << "\n";
					
				}
			}
		}
		
	}else{
		std::cerr << "\nWARNING: No variant signals reported for " << gene << ", with " << n_var << " total variants.\n\n";
		os_log << "\nWARNING: No variant signals reported for " << gene << ", with " << n_var << " total variants.\n\n";
	}
}


void cis_meta_data::conditional_analysis_het(const int& gene_index, std::ostream& os)
{
	
	if( N_CIS[gene_index] <= 0 ){
		std::cerr << "\nWARNING: No variants present for " << gene_id[gene_index] << ".\n";
		std::cerr << "    Note: score.size() == " << score[gene_index].size() << "\n\n";
		
		return;
	}
	
	//cout << "0\n";
	
	std::string& gene = gene_id[gene_index];

	int n_var = N_CIS[gene_index];
	int s_var = S_CIS[gene_index];
	
	Eigen::VectorXd dV( n_var );
	
	// Update IVW weights. Only keep elements from studies where gene is non-missing.
	std::vector<double> ivw_weights;
	
	for( const int& s : study_list[gene_index] ){
		ivw_weights.push_back(ivw[gene_index][s]);
	}
	
	vcov_getter vget(vc, ivw[gene_index], s_var, n_var, study_list[gene_index]);
	
	meta_svar_sumstat meta_ss(vget, ivw_weights);
	
	for( const int& s : study_list[gene_index] ){
		meta_ss.add_sstat(score_perStudy[gene_index][s], dV_perStudy(gene_index, s),  DF_perStudy[gene_index][s], SSR_perStudy[gene_index][s]);
	}
	
	auto snp = [&](const int& i ){ int j = s_var + i; return vc.chr[j] + "_" + std::to_string(vc.pos[j]) + "_" + vc.ref[j] + "_" + vc.alt[j];};
	
	
	// p-value stop threshold is global_opts::LM_ALPHA
	
	std::vector<int> top_snps;
	std::vector<double> acat_stepwise_pvals;
	std::vector<double> svar_stepwise_pvals;
	
	meta_ss.update_meta_ss();
	
	int total_steps = 0;
	while( 1 ){
		
		int n_steps = 0;
		
		int top_snp;
		double svar_stepwise_pval, acat_stepwise_pval;
		
		// std::cout << "BEGIN" << "\t" << gene << "\t";
		// meta_ss.ss_meta.acat_min_pval(top_snp, svar_stepwise_pval, acat_stepwise_pval);
		
		meta_ss.omni_min_pval(top_snp, svar_stepwise_pval, acat_stepwise_pval);
		
		// std::cout << top_snp << "\t" << svar_stepwise_pval << "\t" << acat_stepwise_pval << "\n";
		
		double pval_check = global_opts::step_marginal ? svar_stepwise_pval : acat_stepwise_pval;
		
		// -----------------------------------
		// Forward step. 
		// -----------------------------------
		if( top_snp >= 0 && pval_check >= 0 && ( pval_check < global_opts::LM_ALPHA || top_snps.size() == 0 ) ){
		
			top_snps.push_back(top_snp);
			svar_stepwise_pvals.push_back(svar_stepwise_pval);
			acat_stepwise_pvals.push_back(acat_stepwise_pval);
			
			meta_ss.condition_on_het(top_snp);
			
			n_steps++;
		}
		
		// -----------------------------------
		// Backward step. 
		// -----------------------------------		
		if( global_opts::backward_thresh < 1.00 & top_snps.size() > 1 ){
			
			double max_joint_pvalue = 0;
			int w_rm = -1;
					
			for(int i = 0; i < top_snps.size() - 1; i++){
				double p_hom, p_het, p_aca, p_omn;
				ss_lm_single fm_i = meta_ss.final_model_triform_pval(i, p_hom, p_het, p_aca, p_omn);
				if( p_omn > max_joint_pvalue || p_omn < 0 ){
					max_joint_pvalue = fm_i.pval;
					w_rm = i;
				}
			}
			// update max joint pvalue
			
			if( max_joint_pvalue > global_opts::backward_thresh || max_joint_pvalue < 0 ){
				
				meta_ss.drop_snp(w_rm);
				
				top_snps.erase(top_snps.begin() + w_rm);
				acat_stepwise_pvals.erase(acat_stepwise_pvals.begin() + w_rm);
				svar_stepwise_pvals.erase(svar_stepwise_pvals.begin() + w_rm);
				
				n_steps++;
			}
		}
		
		if( n_steps == 0 ){
			break;
		}else if( top_snps.size() >= global_opts::max_signals ){
			break;
		}
		total_steps++;
		if( total_steps > global_opts::max_steps ){
			std::cerr << "\n\nWARNING: Exceeded max steps. Convergence failed.\n\n";
			break;
		}
	}
	
	
	
	std::string in_studies = std::to_string(study_list[gene_index][0] + 1);
	for( int s = 1; s < study_list[gene_index].size(); s++){
		in_studies += ("," + std::to_string(study_list[gene_index][s] + 1));
	}
	
	// std::cout << "\n\nSTARTING OUTPUT:\n";
	
	for(int i = 0; i < top_snps.size(); i++){
		double p_hom, p_het, p_aca, p_omn;
		double p_hom_m, p_het_m, p_aca_m, p_omn_m;
		ss_lm_single fm_i = meta_ss.final_model_triform_pval(i, p_hom, p_het, p_aca, p_omn);
		ss_lm_single fm_i_m = meta_ss.marginal_triform_pval(i, p_hom_m, p_het_m, p_aca_m, p_omn_m);
		//double marginal_pval_i = meta_ss.ss_meta_0.single_snp_pval(i);
		os << gene << "\t" << in_studies;
		os << "\t" << i + 1 << ":" <<  top_snps.size() << "\t" << snp(top_snps[i]) << "\t" << fm_i.beta << "\t" << fm_i.se << "\t" << p_omn << "\t" << acat_stepwise_pvals[i] << "\t" << p_omn_m << ":" << p_hom_m << "," << p_het_m << "," << p_aca_m << "\t" << svar_stepwise_pvals[i] << "\n";
	}
	// std::cout << "\nEND.\n\n";
}

void cis_meta_data::conditional_analysis(){
	
	// stepwise output

	std::string log_out_name = global_opts::out_prefix + ".cis_meta.stepwise.err_logs";
	std::ofstream os_err(log_out_name.c_str(), std::ofstream::out);

	std::string out_name = global_opts::out_prefix + ".cis_meta.stepwise.tsv";
	std::ofstream os(out_name.c_str(), std::ofstream::out);
	
	std::vector<std::string> col_names{"#gene", "studies", "signal", "variant", "beta", "se", "pval_joint", "pval_signal", "pval_marginal", "pval_stepwise"};
	
	print_header(col_names, os);
	
	// ld buddy output.
	
	std::string buddy_out_name = global_opts::out_prefix + ".cis_meta.buddies.tsv";
	std::ofstream os_b;
	if( global_opts::RSQ_BUDDY < 1 ){
		os_b.open(buddy_out_name.c_str(), std::ofstream::out);
		print_header(std::vector<std::string>{"#signal_variant", "buddy_variant", "r", "rsq"}, os_b);
	}
	
	// loop through genes.
	
	std::cerr << "Completed stepwise meta-analysis for ";
	
	std::string iter_cerr_suffix = " genes out of " + std::to_string(gene_id.size()) + " ...";
	
	print_iter_cerr(1, 0, iter_cerr_suffix);
	for(int i = 0; i < gene_id.size(); ++i){
		int j = i;
		conditional_analysis(j, os, os_b, os_err);
		j = i;
		thinned_iter_cerr(j, i+1, iter_cerr_suffix, 1);
	}
	
	clear_line_cerr();
	std::cerr << "Completed stepwise meta-analysis of " << std::to_string(gene_id.size()) << " total genes.\n";
	
	if( global_opts::RSQ_BUDDY < 1 ) os_b.close();
	os.close();
}


void cis_meta_data::conditional_analysis_het(){

	std::string out_name = global_opts::out_prefix + ".cis_meta.stepwise_het.tsv";
	std::ofstream os(out_name.c_str(), std::ofstream::out);
	
	std::vector<std::string> col_names{"#gene", "studies", "signal", "variant", "beta", "se", "pval_joint", "pval_signal", "pval_marginal", "pval_stepwise"};
	
	print_header(col_names, os);
	
	std::cerr << "Completed heterogeneous stepwise meta-analysis for ";
	
	std::string iter_cerr_suffix = " genes out of " + std::to_string(gene_id.size()) + " ...";
	
	print_iter_cerr(1, 0, iter_cerr_suffix);
	for(int i = 0; i < gene_id.size(); ++i){
		int j = i;
		conditional_analysis_het(j, os);
		j = i;
		thinned_iter_cerr(j, i+1, iter_cerr_suffix, 1);
	}
	
	clear_line_cerr();
	std::cerr << "Completed heterogeneous stepwise meta-analysis of " << std::to_string(gene_id.size()) << " total genes.\n";
	
	os.close();
}


/*int alleles_match(const std::string& ref_0, const std::string& alt_0, const double& freq_0, const std::string& ref_1, const std::string& alt_1, const double& freq_1)
{
	
	//  1 => match
	// -1 => flipped match
	//  0 => no match 
	
	bool is_ambiguous = ambiguous_snv(ref_0, alt_0);

	if( ref_0 == ref_1 && alt_0 == alt_1 ){
		if( std::abs(freq_0 - freq_1) < global_opts::freq_tol )
		{
			return 1;
		}
		if( is_ambiguous )
		{
			if( global_opts::try_match_ambiguous_snv ){
				if( std::abs( (1-freq_0) - freq_1) < global_opts::freq_tol )
				{
					return -1;
				}
			}
		}
	}else if( ref_0 == alt_1 && alt_0 == ref_1 ){
		if( !is_ambiguous )
		{
			if( std::abs( (1-freq_0) - freq_1) < global_opts::freq_tol )
			{
				return -1;
			}
		}else{
			if( global_opts::try_match_ambiguous_snv ){
				//
			}
		}
	}
	
	return 0;
}*/

