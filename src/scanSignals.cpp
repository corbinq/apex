/*  
    Copyright (C) 2020 
    Author: Corbin Quick <qcorbin@hsph.harvard.edu>

    This file is part of YAX.

    YAX is distributed "AS IS" in the hope that it will be 
    useful, but WITHOUT ANY WARRANTY; without even the implied 
    warranty of MERCHANTABILITY, NONINFRINGEMENT, or FITNESS 
    FOR A PARTICULAR PURPOSE.

    The above copyright notice and this permission notice shall 
    be included in all copies or substantial portions of YAX.
*/


#include "Main.hpp"

void scan_signals(bcf_srs_t*& sr, bcf_hdr_t*& hdr,genotype_data& g_data, table& c_data, bed_data& e_data, block_intervals& bm, const bool& rknorm_y, const bool& rknorm_r)
{

	Eigen::MatrixXd &Y = e_data.data_matrix;
	Eigen::MatrixXd &X = c_data.data_matrix;
	
	std::cerr << "Started cis-QTL analysis ...\n";
	
	if( rknorm_y ){
		std::cerr << "Rank-normalizing expression traits ... \n";
		rank_normalize(Y);
	}
	std::cerr << "Scaling expression traits ... \n";
	std::vector<double> y_scale;
	scale_and_center(Y, y_scale);

	Eigen::MatrixXd U = get_half_hat_matrix(X);
	Eigen::VectorXd dV(g_data.var.size());

	if( !global_opts::low_mem ){
		
		std::cerr << "Calculating genotype-covariate covariance...\n";
		
		Eigen::MatrixXd UtG = U.transpose() * g_data.genotypes;
		std::cerr << "Calculating genotype residual variances ...";
		//Eigen::VectorXd SD_vec(UtG.cols());
		for( int i = 0; i < UtG.cols(); ++i)
		{
			
			dV(i) = g_data.genotypes.col(i).squaredNorm() - UtG.col(i).squaredNorm();
			//SD_vec[i] = std::sqrt(g_data.var[i]);
		}
		std::cerr << "Done.\n";
	}
	
	std::cerr << "Calculating expression residuals...\n";
	Eigen::MatrixXd Y_res = resid_from_half_hat(Y, U);
	
	std::cerr << "Scaling expression residuals ...\n";
	scale_and_center(Y_res, e_data.stdev);
	
	if( rknorm_r ){
		std::cerr << "Rank-normalizing expression residuals ...\n";
		rank_normalize(Y_res);
		std::cerr << "Re-residualizing transformed residuals ...\n";
		Eigen::MatrixXd tmp = resid_from_half_hat(Y_res, U);
		Y_res = tmp;
		scale_and_center(Y_res);
	}
	
	// std::cout << Y_res.format(EigenTSV) << "\n";
	// return 0;
	
	Y_res.transposeInPlace();
	
	double n_samples = X.rows();
	double n_covar = X.cols();
	
	// ------------------------------------------
	
	std::string stat_file_path = global_opts::out_prefix + "." + "cis_signal_stats" + ".txt.gz";
	std::string data_file_path = global_opts::out_prefix + "." + "cis_signal_data" + ".txt.gz";
	
	BGZF* stat_file = bgzf_open(stat_file_path.c_str(), "w");
	BGZF* data_file = bgzf_open(data_file_path.c_str(), "w");

	// Write header for stats file 
	std::stringstream os;
	
	std::string stat_header = "#gene\tnsignal\tsnp\tbeta\tse\tp_joint\tp_acat\tp_marginal\tp_sequential\n";
	
	write_to_bgzf(stat_header.c_str(), stat_file);

	// Write header for dataa file 

	os.str(std::string());
	os.clear();
	
	// os << "##RN_Y=" << rknorm_y << ";RN_R=" << rknorm_r << ";NS=" << (int)n_samples << ";NC=" << (int) n_covar << "\n";
	
	os << "#chr" << "\t" << 
		"start" << "\t" << 
		"end" << "\t" << 
		"gene" << "\t" <<
		"snp";
	for( const std::string& sample_id : g_data.ids.keep ){
		os << "\t" << sample_id;
	}
	os << "\n";

	write_to_bgzf(os.str().c_str(), data_file);

	// ------------------------------------------
	// Forward selection and write output
	// ------------------------------------------
	
	int bl = 0;
	
	std::string iter_cerr_suffix = " blocks out of " + std::to_string(bm.size()) + " total";
	std::cerr << "Processed ";
	print_iter_cerr(1, 0, iter_cerr_suffix);
	
	for( int i = 0; i < bm.size(); i++ ){
		
		int n_e = bm.bed_e[i] - bm.bed_s[i] + 1;
		n_e = n_e < Y_res.rows() ? n_e : Y_res.rows()-1;
		
		int n_g = bm.bcf_e[i] - bm.bcf_s[i] + 1;
		n_g = n_g < g_data.n_variants ? n_g : g_data.n_variants-1;
		
		if( n_g > 0 && n_e > 0 && bm.bed_s[i] < Y_res.rows() && bm.bcf_s[i] < g_data.n_variants ){
			
			Eigen::MatrixXd out_mat;
			
			if( global_opts::low_mem ){
				
				g_data.read_genotypes(sr, hdr, bm.bcf_s[i], n_g );
				
				Eigen::SparseMatrix<double>& G = g_data.genotypes;
				
				Eigen::VectorXd UtG_block_sqnm = (U.transpose() * G).colwise().squaredNorm().eval(); // G.middleCols(bm.bcf_s[i], n_g);
				
				for( int si = bm.bcf_s[i], ii = 0; si < bm.bcf_s[i] + n_g; ++si, ++ii)
				{
					dV(si) = G.col(ii).squaredNorm() - UtG_block_sqnm(ii);
					//cout << g_data.var[i] << "\n";
				}
			}
			
			int G_start = bm.bcf_s[i];
			if( global_opts::low_mem ){
				G_start = 0;
			}
			
			// std::cerr << "Get G block\n";
			
			const Eigen::SparseMatrix<double>& G = g_data.genotypes.middleCols(G_start, n_g);
			// std::cerr << g_data.genotypes.rows() << "\t" << g_data.genotypes.cols() << "\n";
			
			// std::cerr << "Get outmat\n";
			
			out_mat = (Y_res.middleRows(bm.bed_s[i], n_e)* G).transpose().eval();
			
			// std::cerr << "Get UtG\n";
			
			Eigen::MatrixXd UtG = U.transpose() * G;
			
			for( int jj = bm.bed_s[i], jm = 0; jj < bm.bed_s[i] + n_e; jj++, jm++){
				
				// std::cerr << "Start U V\n";
				
				const Eigen::VectorXd& U = out_mat.col(jm);
				const Eigen::VectorXd& V = dV.segment(bm.bcf_s[i], n_g);
				
				// std::cerr << "out_mat: " << out_mat.rows() << "\t" << out_mat.cols() << "\n";
				
				
				// std::cerr << U.size() << "\t" << V.size() << "\n";
				
				// std::cerr << "Start vget\n";
			
				indiv_vcov_getter vget(G, UtG);
				
				// std::cerr << G.rows() << "\t" << G.cols() << "\n";
				// std::cerr << UtG.rows() << "\t" << UtG.cols() << "\n";
				
				double y_sc = y_scale[jj] * e_data.stdev[jj];
				
				// std::cerr << "Start forward\n";
				
				double SSR = e_data.stdev[jj] * e_data.stdev[jj] * (n_samples - 1);
				double IVW =  (n_samples - n_covar)/SSR;
				// double adj = e_data.stdev[jj] * IVW;
				
				double adj = std::sqrt(n_samples - n_covar)/std::sqrt(n_samples - 1);
				
				forward_lm flm(adj * U, V, n_samples, n_samples - n_covar, e_data.stdev[jj]/adj, vget, global_opts::LM_ALPHA);
				
				// std::cerr << "End forward\n";
				
				int bcf_s_b = bm.bcf_s[i];
		
				auto snp = [&](const int& kk ){ int j = bcf_s_b + kk; return clean_chrom(g_data.chr[j]) + "_" + std::to_string(g_data.pos[j]) + "_" + g_data.ref[j] + "_" + g_data.alt[j];};
				
				if( flm.beta.size() > 0 ){
					
					std::stringstream os;
					
					for(int bi = 0; bi < flm.beta.size(); bi++)
					{
						double slope = flm.beta[bi];
						if( g_data.flipped[bcf_s_b + flm.keep[bi]] ){
							slope *= -1.00; 
						}
						// os.precision(4);
						os << e_data.gene_id[jj];
						os << "\t" << bi+1 << ":" << flm.beta.size() << "\t" << snp(flm.keep[bi]) << "\t" << 
						slope << "\t" << flm.se[bi] << "\t" << flm.pval_joint[bi]  << "\t" << flm.pval_adj[bi] << "\t";
						// os.precision(2);
						os << flm.pval_0[bi] << "\t" << flm.pval_seq[bi] << "\n";
					}
					
					write_to_bgzf(os.str().c_str(), stat_file);
					
					os.str(std::string());
					os.clear();

					for(const int& k : flm.keep ){
						os << clean_chrom(e_data.chr[jj]) << "\t" << 
							e_data.start[jj] << "\t" << 
							e_data.end[jj] << "\t" << 
							e_data.gene_id[jj] << "\t" << 
							snp(k);
						Eigen::VectorXd G_k = G.col(k);
						for(int bi = 0; bi < G_k.size(); bi++){
							os << "\t" << G_k(bi);
						}
						os << "\n";
					}
					
				
					write_to_bgzf(os.str().c_str(), data_file);
				}
			}
			
		}else{
			std::cerr << "\nERROR: " << bl << "; " << bm.bed_s[i] << ", " << n_e << "; " << bm.bcf_s[i] << ", " << n_g << "\n"; 
			abort();
		}
		
		print_iter_cerr(bl, bl+1, iter_cerr_suffix);
		
		bl++;
	}
	std::cerr << "\n";

	bgzf_close(stat_file);
	bgzf_close(data_file);

	return;
}
