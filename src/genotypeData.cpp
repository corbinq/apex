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

#include "genotypeData.hpp"

bool sp_geno_fmt = true;

void read_sparse_GRM(const std::string& filename, Eigen::SparseMatrix<double>& GRM, const std::vector<std::string>& kp_ids, const double& r_scale, const int& r_col, std::vector<int>& related)
{
	int n = kp_ids.size();
	
	GRM = Eigen::SparseMatrix<double>(n,n);
	
	if( filename == "" ){
		GRM.setIdentity();
		return;
	}
	
	std::unordered_map<std::string, int> id_map;
	for( int i = 0; i < n; i++ ){
		id_map[ kp_ids[i] ] = i;
	}
	
	std::vector<std::string> id1, id2;
	std::vector<double> val;
	
	data_parser dp;
	dp.add_field(id1, 0);
	dp.add_field(id2, 1);
	dp.add_field(val, r_col - 1);
	
	int nr = 0;
	dp.parse_file(filename, nr);

	using td = Eigen::Triplet<double>;

	std::vector<td> triplets;
	
	bool add_diag = true;
	
	for( int i = 0; i < nr; i++ )
	{
		const auto f1 = id_map.find(id1[i]);
		const auto f2 = id_map.find(id2[i]);
		
		if( add_diag ){
			if( id1[i] == id2[i] ){
				add_diag = false;
			}
		}
		
		if( f1 != id_map.end() && f2 != id_map.end() )
		{
			int& ii = f1->second;
			int& jj = f2->second;
			if( ii != jj ){
				related.push_back(ii);
				related.push_back(jj);
			}
			triplets.push_back(td(ii,jj,r_scale*val[i]));
			triplets.push_back(td(jj,ii,r_scale*val[i]));
		}
		
	}
	if( add_diag ){
		for(int i = 0; i < n; ++i) triplets.push_back(td(i,i,1.0));
	}
	
	std::sort( related.begin(), related.end() );
	related.erase( std::unique( related.begin(), related.end() ), related.end() );
	
	GRM.setFromTriplets(triplets.begin(), triplets.end());
}

void genotype_data::read_bcf_variants(bcf_srs_t*& sr, bcf_hdr_t*& hdr, int& n_var, bool store_geno, bool scan_geno)
{

	if( store_geno ){
		std::cerr << "Processed genotype data for ";
	}else{
		std::cerr << "Processed variant data for ";
	}
	std::string iter_cerr_suffix = " variants ... ";
	
	print_iter_cerr(1, 0, iter_cerr_suffix);
	int last = 0;
	
	initialize_genotypes(n_var);
	n_var = 0;
	geno_size = 0;
	
	while( bcf_sr_next_line(sr) )
	{
		bcf1_t *rec = bcf_sr_get_line(sr,0);	
		if( process_bcf_variant(rec, hdr, store_geno, scan_geno) ){
			n_var++;
			if( store_geno ) geno_size++;
		}
		thinned_iter_cerr(last, n_var, iter_cerr_suffix, 2500);
	}
	print_iter_cerr(last, n_var, iter_cerr_suffix);
	
	geno_start = 0;

	get_ld_index();
	
	return;
}

int genotype_data::get_ld_index()
{
	int bp = global_opts::ld_window_bp;
	
	int max_entries = 0;
	
	int ii = 0;
	int ii_s = 0;
	int ii_e = 0;
	while( ii < pos.size() ){
		while( pos[ii] - pos[ii_s] > bp && ii_s < pos.size() ){
			ii_s++;
		}
		while( pos[ii_e] - pos[ii] <= bp && ii_e < pos.size() ){
			ii_e++;
		}
		
		ld_index_range.push_back( std::make_pair(ii, ii_e) );
		
		if( ii_e - ii - 1 > max_entries ){
			max_entries = ii_e - ii - 1;
		}
		
		ii++;
	}
	return max_entries;
}

void genotype_data::clear_genotypes()
{
	genotypes.setZero();
	genotypes.data().squeeze(); 
}

void genotype_data::clear()
{
	n_variants = 0;
	
	chr.clear();
	pos.clear();
	ref.clear();
	alt.clear();
	
	mean.clear();
	var.clear();
	
	flipped.clear();
	
	clear_genotypes();
}

inline bool genotype_data::add_bcf_genotypes(int*& gt_rec, const int& col_n, double& mean_, double& var_, bool& flip_, const bool store_geno)
{
	int n = 0;
	mean_ = 0;
	var_ = 0;
	flip_ = false;
	sparse_gt gts;
	std::vector<int> missing;
	if( sp_geno_fmt ){
		gts.resize(ids.idx.size());
	}else{
		if( dense_genotypes.rows() < ids.keep.size() ){
			dense_genotypes.conservativeResize(ids.keep.size(), Eigen::NoChange);
		}
		if( dense_genotypes.cols() < col_n ){
			dense_genotypes.conservativeResize(Eigen::NoChange, col_n + chunk_size);
		}
	}
	int n_obs = 0;
	
	int n_2 = 0;
	int n_1 = 0;
	
	for(const int& i: ids.idx){

		int a1 = bcf_gt_allele(gt_rec[i*2+0]);
		int a2 = bcf_gt_allele(gt_rec[i*2+1]);
		if(a1 < 0 || a2 < 0)
		{
			if( global_opts::exclude_missing > 0 ){
				return false;
			}else{
				// a1 = a1 < 0 ? 0 : a1;
				// a2 = a2 < 0 ? 0 : a2;
				if(store_geno){
					if( sp_geno_fmt ){
						gts.set_gt(n,-1);
					}else{
						missing.push_back(i);
					}
				}
			}
		}else{
		
			if( a1 + a2 > 0 ){
				a1 += a2;
				if( a1 > 1){
					n_2++;
				}else{
					n_1++;
				}
				if(store_geno){
					if( sp_geno_fmt ){
						gts.set_gt(n,a1);
					}else{
						// DO DENSE
						// std::cerr << "DO DENSE";
						dense_genotypes(n, col_n) = a1;
					}
				}
				mean_ += a1;
				var_ += a1*a1;
			}
			n_obs++;
		}
		n++;
	}
	
	mean_ = mean_/((double) n_obs);
	var_ =  (var_ - n_obs*mean_*mean_)/ ( (double) n_obs + 1 );
	
	if(  n_2 > n_obs - n_1 - n_2  ){
		flip_ = true;
		mean_ = 2.0 - mean_;
		if( store_geno && sp_geno_fmt ) gts.flip(n);
	}
	if( store_geno && sp_geno_fmt ){
		// checkResize(col_n);
		gts.add_gt_sparsemat(genotypes, col_n);
	}
	if( store_geno && !sp_geno_fmt && missing.size() > 0 ){
		for(const int& i : missing ){
			dense_genotypes(i, col_n) = mean_;
		}
	}
	
	return true;
}

inline bool genotype_data::add_bcf_dosages(float*& ds_rec, const int& col_n, double& mean_, double& var_, bool& flip_, const bool store_geno)
{
	int n = 0;
	int n_obs = 0;
	mean_ = 0;
	var_ = 0;
	flip_ = false;
	sparse_ds sp_ds;
	std::vector<int> missing;
	if( sp_geno_fmt ){
		sp_ds.resize(ids.idx.size());
	}else{
		if( dense_genotypes.rows() < ids.keep.size() ){
			dense_genotypes.conservativeResize(ids.keep.size(), Eigen::NoChange);
		}
		if( dense_genotypes.cols() < col_n ){
			dense_genotypes.conservativeResize(Eigen::NoChange, col_n + chunk_size);
		}
	}
	
	int n_2 = 0;
	int n_0 = 0;
	
	for(const int& i: ids.idx)
	{
		float ds = ds_rec[i];
		if( ds < 0 ){
			if(store_geno){
				if( sp_geno_fmt ){
					sp_ds.set_ds(n,-1.00);
				}else{
					missing.push_back(i);
				}
			}
		}else{
			if( ds <= global_opts::dosage_thresh ){
				ds = 0.00;
				n_0++;
			}else if( ds >= 2.00 - global_opts::dosage_thresh ){
				ds = 2.00;
				n_2++;
			}
			if(store_geno){
				if( sp_geno_fmt ){
					sp_ds.set_ds(n,ds);
				}else{
					// std::cerr << "DO DENSE";
					dense_genotypes(n, col_n) = ds;
				}
			} 
			mean_ += ds;
			var_ += ds*ds;
			n_obs++;
		}
		n++;
	}
	
	mean_ = mean_/((double) n_obs);
	var_ =  (var_ - n_obs*mean_*mean_)/ ( (double) n_obs + 1 );
	
	if( n_2 > n_0 ){
		flip_ = true;
		mean_ = 2.0 - mean_;
		if( store_geno && sp_geno_fmt ) sp_ds.flip();
	}
	if( store_geno && sp_geno_fmt ){
		// checkResize(col_n);
		sp_ds.add_ds_sparsemat(genotypes, col_n);
	}
	if( store_geno && !sp_geno_fmt && missing.size() > 0 ){
		for(const int& i : missing ){
			dense_genotypes(i, col_n) = mean_;
		}
	}
	
	return true;
}

void genotype_data::read_bcf_header(bcf_hdr_t* hdr){

	n_samples = bcf_hdr_nsamples(hdr);
	
	for(int i = 0; i < n_samples; i++){
		ids.file.push_back(hdr->samples[i]);
	}
	if( ids.keep.size() == 0 ){
		ids.setKeepIDs(ids.file);
	}
}

void genotype_data::freeze_genotypes(){
	genotypes.conservativeResize(ids.keep.size(), geno_size);
	genotypes.finalize();
	genotypes.makeCompressed();
}

inline bool genotype_data::process_bcf_variant(bcf1_t*& rec, bcf_hdr_t*& hdr, bool store_geno, bool scan_geno){
	
	if( rec->n_allele != 2 ){
		std::cerr << "failed rec->n_allele" << rec->n_allele << "\n";
		return false;
	}
	
	if( n_samples <= 0 ){
		if( ids.keep.size() == 0 ){
			read_bcf_header(hdr);
		}else{
			n_samples = ids.keep.size();
		}
	}
	
	double mean_ = -1; 
	double var_ = -1;
	
	bool flip_ = false;
	bool keep_ = true;
	
	// std::string snp_id(rec->d.id);
	
	if( scan_geno ){
		
		if( global_opts::use_dosages ){
			int nds_arr = 0;
			float *ds = NULL;
			bcf_get_format_float(hdr, rec, "DS", &ds, &nds_arr)/n_samples;
			keep_ = add_bcf_dosages(ds, n_variants, mean_, var_, flip_, store_geno);
			delete ds;
		}else{
			int ngt_arr = 0;
			int *gt = NULL;
			bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt_arr)/n_samples;
			keep_ = add_bcf_genotypes(gt, n_variants, mean_, var_, flip_, store_geno);
			delete gt;
		}
		
		if( !keep_ ){
			//genotypes.conservativeResize(n_variants, Eigen::NoChange);
			std::cout << "FAILED\n";
			return false;
		}
		
	}
	
	n_variants++;
		
	chr.push_back(bcf_hdr_id2name(hdr, rec->rid));
	pos.push_back(rec->pos + 1);
	//rsid.push_back(rsid_);
	ref.push_back(rec->d.allele[0]);
	alt.push_back(rec->d.allele[1]);
	
	flipped.push_back(flip_);
	keep.push_back(keep_);
	
	mean.push_back(mean_);
	var.push_back(var_);
	
	return true;
}

bool genotype_data::record_matches(bcf1_t*& rec, bcf_hdr_t*& hdr, const int& i){
	//cerr <<  bcf_hdr_id2name(hdr, rec->rid) << "\t" << (rec->pos + 1) << "\t" 
	if( bcf_hdr_id2name(hdr, rec->rid) == chr[i] ){
		if( (rec->pos + 1) == pos[i] ){
			if( rec->d.allele[0] == ref[i] ){
				if( rec->d.allele[1] == alt[i] ){
					return true;
				}
			}
		}
	}
	return false;
}

void genotype_data::read_genotypes(bcf_srs_t*& sr, bcf_hdr_t*& hdr, const int& i_s, const int& n_s)
{
	clear_genotypes();
	resize_genotypes(n_s);
	
	int i_e = i_s + n_s - 1;
	if( chr[i_s] != chr[i_e] ){
		std::cerr << "Fatal: regional chr mismatch in genotype_data::read_genotypes. \n";
		exit(1);
	}
	bcf_seek(sr, chr[i_s], pos[i_s]);
	
	int i_i = i_s;
	
	geno_start = i_s;
	geno_size = 0;

	int r_i = 0;
	
	while( bcf_sr_next_line(sr) && r_i < n_s )
	{
		bcf1_t *rec = bcf_sr_get_line(sr,0);
		
		if( record_matches(rec, hdr, i_i) ){
			
			if( keep[i_i] ){

				bool is_flipped = false;
				
						
				if( global_opts::use_dosages ){
					int nds_arr = 0;
					float *ds = NULL;
					bcf_get_format_float(hdr, rec, "DS", &ds, &nds_arr)/n_samples;
					keep[i_i] = add_bcf_dosages(ds, r_i, mean[i_i], var[i_i], is_flipped, true);
					delete ds;
				}else{
					int ngt_arr = 0;
					int *gt = NULL;
					int n_gts = bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt_arr)/n_samples;
					keep[i_i] = add_bcf_genotypes(gt, r_i, mean[i_i], var[i_i], is_flipped, true);
					delete gt;
				}
					
				flipped[i_i] = is_flipped;
				
				
			}

			r_i++;
			i_i++;
			geno_size++;
		}
	}
	if( i_i < i_e ){
		std::cerr << "Fatal: Failed to fill region in genotype_data::read_genotypes. \n";
		std::cerr << i_s << ", " << i_i << ", " << i_e << "\n";
		exit(1);
	}
	freeze_genotypes();
}
