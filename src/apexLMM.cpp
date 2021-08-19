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


#include "apexLMM.hpp"

void fit_LMM_null_models(table& c_data, bed_data& e_data, Eigen::SparseMatrix<double>& L, Eigen::VectorXd& GRM_lambda, const bool& rknorm_y, const bool& rknorm_r)
{
	
	Eigen::SparseMatrix<double> Q;
	Eigen::VectorXd Q_lambda;

	subset_eigen(L, GRM_lambda, Q, Q_lambda);
	
	std::cerr << "Reordering trait and covariate matrices ...\n";
	
	Eigen::MatrixXd& Y = e_data.data_matrix;
	Eigen::MatrixXd& C = c_data.data_matrix;
	
	double n_traits = Y.cols();
	double n_samples = Y.rows();
	double n_covar = C.cols();

	// std::cerr << "Started cis-QTL analysis ...\n";
	
	if( rknorm_y ){
		std::cerr << "Rank-normalizing expression traits ... \n";
		rank_normalize(Y);
	}
	std::cerr << "Scaling expression traits ... \n";
	std::vector<double> y_scale;
	scale_and_center(Y, y_scale);
	
	std::cerr << "Calculating partial rotations ...\n";
	Eigen::MatrixXd QtC = (Q.transpose() * C).eval();
	Eigen::MatrixXd QtY = (Q.transpose() * Y).eval();
	Eigen::MatrixXd CtY = (C.transpose() * Y).eval();
	Eigen::MatrixXd CtC = (C.transpose() * C).eval();
	Eigen::MatrixXd CtC_i = (CtC.inverse()).eval();
	
	std::cerr << "Rotating expression and covariates ... ";
	//Y = (L.transpose() * Y).eval();
	Eigen::MatrixXd X = (L.transpose() * C).eval();
	std::cerr << "Done.\n";
	
	std::string theta_file_path = global_opts::out_prefix + "." + "theta" + ".gz";
	
	BGZF* theta_file;
	theta_file = bgzf_open(theta_file_path.c_str(), "w");
	
	std::string iter_cerr_suffix = " traits out of " + std::to_string(Y.cols()) + " total";
	std::cerr << "Fit null model for ";
	print_iter_cerr(1, 0, iter_cerr_suffix);
	
	std::vector<double> phi_v(Y.cols());
	std::vector<double> hsq_v(Y.cols());
	std::vector<double> sigma_v(Y.cols());
	std::vector<double> SSR_v(Y.cols());
	
	e_data.stdev.resize(Y.cols());
	
	for(int j = 0; j < Y.cols(); j++ ){

		Y.col(j) = (L.transpose() * Y.col(j)).eval();

		LMM_fitter fit(X, Y.col(j), GRM_lambda);
		fit.fit_REML();
		
		const DiagonalXd& Vi = fit.Vi;
		const double& sigma2 = fit.sigma2;
		const double& phi = fit.phi;
		
		
		
		// DiagonalXd Vi = Eigen::VectorXd::Ones(Y.col(j).size()).asDiagonal();
		// double sigma2 = 1.00;
		// double phi = 0.05;
		
		double tau2 = phi*sigma2;
		double hsq = tau2 / (tau2 + sigma2);
		double scale = std::sqrt(tau2 + sigma2);
		
		// Calculate "sum of squared residuals"

		DiagonalXd Psi = calc_Psi(phi, Q_lambda);
		Eigen::MatrixXd XtDX = (CtC - QtC.transpose() * Psi * QtC )/(1.00 + phi);
		Eigen::VectorXd XtDy = X.transpose() * Vi * Y.col(j);
		Eigen::VectorXd b = XtDX.colPivHouseholderQr().solve(XtDy);
		Eigen::VectorXd y_hat = X * b;
		
		Eigen::VectorXd y_res = Y.col(j) - y_hat;
		
		SSR_v[j] = y_res.dot(Vi * Y.col(j))/sigma2;
	
		
		// if( !e_data.is_residual ){
			// DiagonalXd Psi = calc_Psi(phi, Q_lambda);
			// Eigen::MatrixXd XtDX = (CtC - QtC.transpose() * Psi * QtC )/(1.00 + phi);
			// Eigen::VectorXd XtDy = X.transpose() * Vi * Y.col(j);
			// Eigen::VectorXd b = XtDX.colPivHouseholderQr().solve(XtDy);
			// Eigen::VectorXd y_hat = X * b;
			
			// Eigen::VectorXd y_res = Y.col(j) - y_hat;
			
			// SSR_v[j] = y_res.dot(Vi * Y.col(j))/sigma2;
			
			// // Y now stores the rotated residuals. 
			// Y.col(j) = (Vi * y_res/std::sqrt(sigma2)).eval();
		// }
		
		std::stringstream theta_line;
		
		theta_line <<
			clean_chrom(e_data.chr[j]) << "\t" << 
			e_data.start[j] << "\t" << 
			e_data.end[j] << "\t" << 
			e_data.gene_id[j] << "\t" << 
			sigma2 << "\t" << 
			tau2 << "\t" << 
			phi << "\t" <<
			SSR_v[j] << "\t" <<
			y_scale[j] << "\n";
		
		write_to_bgzf(theta_line.str().c_str(), theta_file);
			
		if(  global_opts::write_resid_mat ){
		
			Y.col(j) = (Vi * y_res/std::sqrt(sigma2)).eval();
			
			Y.col(j) = (L * Y.col(j)).eval();
		}
		
		
		print_iter_cerr(j, j+1, iter_cerr_suffix);
	}
	
	
	bgzf_close(theta_file);
	build_tabix_index(theta_file_path, 1);
	
	
	if(  global_opts::write_resid_mat ){
		// Y = (L * Y).eval();
		std::string bed_out = global_opts::out_prefix + ".lmm_resid.bed.gz";
		std::string bed_header = "## APEX_LMM_RESID";
		e_data.write_bed(bed_out, bed_header);
	}
	
}


void fit_LMM_null_models_low_rank(const int& n_fac, table& c_data, bed_data& e_data, const bool& rknorm_y, const bool& rknorm_r, Eigen::MatrixXd& Y_epc )
{
	
	Eigen::MatrixXd Q;
	Eigen::VectorXd Q_lambda;
	
	Eigen::MatrixXd& Y = e_data.data_matrix;
	Eigen::MatrixXd& C = c_data.data_matrix;
		
	double n_traits = Y.cols();
	double n_samples = Y.rows();
	double n_covar = C.cols();

	if( rknorm_y ){
		std::cerr << "Rank-normalizing expression traits ... \n";
		rank_normalize(Y);
		rank_normalize(Y_epc);
	}
	std::cerr << "Scaling expression traits ... \n";
	std::vector<double> y_scale;
	scale_and_center(Y, y_scale);
	scale_and_center(Y_epc);
	
	// bool resid_ePC = true;
	
	// if( resid_ePC ){
			
		// Eigen::MatrixXd U = get_half_hat_matrix(C);

		// std::cerr << "Calculating expression residuals...\n";
		// make_resid_from_half_hat(Y_epc, U);
	// }
	
	std::cerr << "Estimating latent factors ... \n";
	calc_eGRM_PCs(Q, Q_lambda, Y_epc, n_fac);
	std::cerr << "Done.\n";
	
	Y_epc.resize(0,0);
	
	std::cerr << "Calculating partial rotations ...\n";
	Eigen::MatrixXd QtC = (Q.transpose() * C).eval();
	Eigen::MatrixXd QtY = (Q.transpose() * Y).eval();
	Eigen::MatrixXd CtY = (C.transpose() * Y).eval();
	Eigen::MatrixXd CtC = (C.transpose() * C).eval();
	Eigen::MatrixXd CtC_i = (CtC.inverse()).eval();
	
	std::vector<double> hsq_vals{0.0, 0.5, 1.0};
	
	std::string theta_file_path = global_opts::out_prefix + "." + "theta" + ".gz";
	BGZF* theta_file;
	theta_file = bgzf_open(theta_file_path.c_str(), "w");
	
	// -------------------------------------
	// Begin fitting null models. 
	// -------------------------------------
	
	std::vector<double> phi_v(Y.cols());
	std::vector<double> hsq_v(Y.cols());
	std::vector<double> sigma_v(Y.cols());
	std::vector<double> SSR_v(Y.cols());
	
	// Eigen::MatrixXd PY( Y.rows(), Y.cols() );
	
	e_data.stdev.resize(Y.cols());
	
	theta_data t_data;
	bool use_theta = false;
	
	// if( theta_path != "" ){
		// std::cerr << "Set null model for ";
		// t_data.open(theta_path);
		// use_theta = true;
	// }else{
		std::cerr << "Fit null model for ";
	// }
	
	std::string iter_cerr_suffix = " traits out of " + std::to_string(Y.cols()) + " total";
	print_iter_cerr(1, 0, iter_cerr_suffix);
	
	int last_j = 0;
	for(int j = 0; j < Y.cols(); j++ ){
		DiagonalXd Vi;
		double sigma2, phi, tau2;
		
		double yty = Y.col(j).squaredNorm();
		double sample_size = Y.rows();
		
		LMM_fitter_low_rank fit(CtC, QtC, CtY.col(j), QtY.col(j), yty, Q_lambda, sample_size);
		
		fit.fit_REML();
		
		// Vi = fit.Vi;
		sigma2 = fit.sigma2;
		phi = fit.phi;
		tau2 = phi*sigma2;
		
		std::stringstream theta_line;
		
		theta_line <<
			clean_chrom(e_data.chr[j]) << "\t" << 
			e_data.start[j] << "\t" << 
			e_data.end[j] << "\t" << 
			e_data.gene_id[j] << "\t" << 
			sigma2 << "\t" << 
			tau2 << "\t" << 
			phi << "\n";
		
		write_to_bgzf(theta_line.str().c_str(), theta_file);
		
		if ( global_opts::write_resid_mat ){
			
			double hsq = tau2 / (tau2 + sigma2);
			double scale = std::sqrt(tau2 + sigma2);
			
			DiagonalXd Psi = calc_Psi_low_rank(phi, Q_lambda);
			Eigen::MatrixXd XtDX = (CtC - QtC.transpose() * Psi * QtC );
			
			Eigen::VectorXd XtDy = QtC.transpose() * Psi * QtY.col(j);
			XtDy = (CtY.col(j) - XtDy).eval();
			Eigen::VectorXd b = XtDX.colPivHouseholderQr().solve(XtDy);
			Eigen::VectorXd y_fixeff = C * b;
			
			Eigen::VectorXd v1 = QtC * b;
			Eigen::VectorXd v2 = QtY.col(j) - v1;
			Eigen::VectorXd y_raneff = Q * Psi * v2;
			Eigen::VectorXd y_resid = Y.col(j) - y_raneff - y_fixeff;
			
			SSR_v[j] = (yty - QtY.col(j).dot(Psi * QtY.col(j)) -  XtDy.dot(b))/sigma2;
			
			Y.col(j) = (y_resid/std::sqrt(sigma2)).eval();
			
			phi_v[j] = phi;
			hsq_v[j] = hsq;
			sigma_v[j] = sigma2;
			
			e_data.stdev[j] = std::sqrt(sigma2);
			
		}
		print_iter_cerr(j, j+1, iter_cerr_suffix);
	}
	std::cerr << "\n";
	
	bgzf_close(theta_file);
	build_tabix_index(theta_file_path, 1);
	
	if(  global_opts::write_resid_mat ){
		std::string bed_out = global_opts::out_prefix + ".lmm_resid.bed.gz";
		std::string bed_header = "## APEX_LMM_RESID";
		e_data.write_bed(bed_out);
	}
}



void calculate_V_anchor_points( Eigen::MatrixXd& V_mat, genotype_data& g_data, const Eigen::MatrixXd& C, const std::vector<double>& hsq_vals, const Eigen::MatrixXd& CtC, const Eigen::MatrixXd& CtC_i, const Eigen::MatrixXd& QtG, Eigen::MatrixXd& QtC, const Eigen::VectorXd Q_lambda ){
	
	BGZF* anchor_file;
	
	int n_hsq = hsq_vals.size();
	int n_snps = g_data.genotypes.cols();
	int n_samples = g_data.genotypes.rows();
	int n_cov = CtC.rows();
	
	if( global_opts::write_v_anchors ){
		std::string anchor_path = global_opts::out_prefix + ".gvar_points.gz";
		anchor_file = bgzf_open(anchor_path.c_str(), "w");
		
		write_to_bgzf("##\t" + std::to_string(n_samples) + "\t" + std::to_string(n_snps) + "\t" + std::to_string(n_cov) + "\n", anchor_file);
		
		write_to_bgzf("#chrom\tpos\tref\talt", anchor_file);
		
		std::stringstream out_line;
		for(const double& val : hsq_vals){
			out_line << "\tVR" << 100.0*val;
		}
		out_line << "\n";
		write_to_bgzf(out_line.str().c_str(), anchor_file);
	}
	
	std::vector<Eigen::MatrixXd> M_list;
	V_mat = Eigen::MatrixXd(n_snps, n_hsq);
	
	for( const double& hsq : hsq_vals ){
		M_list.push_back(
			(CtC - (QtC.transpose() * calc_Psi(hsq, Q_lambda) * QtC)).inverse()
		);
	}
	
	std::cerr << "Calculating genotype-covariate covariance ... \n";
	Eigen::MatrixXd CtG = (C.transpose() * g_data.genotypes).eval();
	
	std::vector<double> GtG_diag(n_snps);		
	
	std::cerr << "Calculating genotype residual variances ...";

	for( int i = 0; i < n_snps; i++)
	{
		GtG_diag[i] = g_data.genotypes.col(i).squaredNorm();
		g_data.var[i] = GtG_diag[i] - CtG.col(i).dot(CtC_i * CtG.col(i));
		V_mat(i,0) = g_data.var[i];
	}
	std::cerr << "Done.\n";
	
	for(int j = 1; j < n_hsq; j++){
		
		std::cerr << "Calculating genotype variance anchor points (" << j+1 << "/" << n_hsq << ") ... \n";

		Eigen::MatrixXd Psi_a = calc_Psi(hsq_vals[j], Q_lambda) * QtG;
		Eigen::MatrixXd Res = CtG - QtC.transpose() * Psi_a;
		
		for( int i = 0; i < n_snps; i++)
		{
			V_mat(i,j) = (
				GtG_diag[i] - QtG.col(i).dot(Psi_a.col(i)) - Res.col(i).dot(M_list[j] * Res.col(i))
			);
		}
		
	}
	
	if( global_opts::write_v_anchors ){
		std::cerr << "Writing genotype variance anchor points ... ";
		for( int i = 0; i < n_snps; i++)
		{
			std::stringstream out_line;
			out_line << 
				clean_chrom(g_data.chr[i]) << "\t" <<
				g_data.pos[i] << "\t" <<
				g_data.ref[i] << "\t" <<
				g_data.alt[i];
			for(int j = 0; j < n_hsq; j++){
				out_line << "\t" << V_mat(i,j);
			}
			out_line << "\n";
			write_to_bgzf(out_line.str().c_str(), anchor_file);
		}
		bgzf_close(anchor_file);
	}
	
	V_mat = (V_mat * (getPredParams(hsq_vals).transpose())).eval();
	
	std::cerr << "Done.\n";
	
	return;
}


void calculate_V_anchor_points_low_rank( Eigen::MatrixXd& V_mat, genotype_data& g_data, const Eigen::MatrixXd& C, const std::vector<double>& hsq_vals, const Eigen::MatrixXd& CtC, const Eigen::MatrixXd& CtC_i, const Eigen::MatrixXd& QtG, Eigen::MatrixXd& QtC, const Eigen::VectorXd Q_lambda ){
	
	
	BGZF* anchor_file;
	std::string anchor_path = global_opts::out_prefix + "." + "gvar_points.gz";
	
	int n_hsq = hsq_vals.size();
	int n_snps = g_data.genotypes.cols();
	int n_samples = g_data.genotypes.rows();
	int n_cov = CtC.rows();
	
	if( global_opts::write_v_anchors ){
		
		anchor_file = bgzf_open(anchor_path.c_str(), "w");
		
		write_to_bgzf("##\t" + std::to_string(n_samples) + "\t" + std::to_string(n_snps) + "\t" + std::to_string(n_cov) + "\n", anchor_file);
		
		write_to_bgzf("#chrom\tpos\tref\talt", anchor_file);
		std::stringstream out_line;
		for(const double& val : hsq_vals){
			out_line << "\tVR" << 100.0*val;
		}
		out_line << "\n";
		write_to_bgzf(out_line.str().c_str(), anchor_file);
	}
	
	std::vector<Eigen::MatrixXd> M_list;
	V_mat = Eigen::MatrixXd(n_snps, n_hsq);
	
	for( const double& hsq : hsq_vals ){
		M_list.push_back(
			(CtC - (QtC.transpose() * calc_Psi_low_rank(hsq, Q_lambda) * QtC)).inverse()
		);
	}
	
	std::cerr << "Calculating genotype-covariate covariance ... \n";
	Eigen::MatrixXd CtG = (C.transpose() * g_data.genotypes).eval();
	
	std::vector<double> GtG_diag(n_snps);		
	
	std::cerr << "Calculating genotype residual variances ...";

	for( int i = 0; i < n_snps; i++)
	{
		GtG_diag[i] = g_data.genotypes.col(i).squaredNorm();
		g_data.var[i] = GtG_diag[i] - CtG.col(i).dot(CtC_i * CtG.col(i));
		V_mat(i,0) = g_data.var[i];
	}
	std::cerr << "Done.\n";
	
	for(int j = 1; j < n_hsq; j++){
		std::cerr << "Calculating genotype variance anchor points (" << j+1 << "/" << n_hsq << ") ... \n";

		Eigen::MatrixXd Psi_a = calc_Psi_low_rank(hsq_vals[j], Q_lambda) * QtG;
		Eigen::MatrixXd Res = CtG - QtC.transpose() * Psi_a;
		
		for( int i = 0; i < n_snps; i++)
		{
			V_mat(i,j) = (
				GtG_diag[i] - QtG.col(i).dot(Psi_a.col(i)) - Res.col(i).dot(M_list[j] * Res.col(i))
			);
		}
	}
	
	if( global_opts::write_v_anchors ){
		std::cerr << "Writing genotype variance anchor points ... ";
		for( int i = 0; i < n_snps; i++)
		{
			std::stringstream out_line;
			out_line << 
				clean_chrom(g_data.chr[i]) << "\t" <<
				g_data.pos[i] << "\t" <<
				g_data.ref[i] << "\t" <<
				g_data.alt[i];
			for(int j = 0; j < n_hsq; j++){
				out_line << "\t" << V_mat(i,j);
			}
			out_line << "\n";
			write_to_bgzf(out_line.str().c_str(), anchor_file);
		}
		bgzf_close(anchor_file);
		build_tabix_index(anchor_path, 0);
	}
	
	V_mat = (V_mat * (getPredParams(hsq_vals).transpose())).eval();
	
	std::cerr << "Done.\n";
	
	return;
}

void read_V_anchor_points( const std::string& anchor_path, Eigen::MatrixXd& V_mat, genotype_data& g_data, const std::vector<double>& hsq_vals)
{
	
	int n_hsq = hsq_vals.size();
	int n_snps = g_data.genotypes.cols();
	
	V_mat = Eigen::MatrixXd(n_snps, n_hsq);
	
	std::vector<std::string> chr;
	std::vector<int> pos;
	std::vector<std::string> ref;
	std::vector<std::string> alt;
	
	data_parser dp;
	
	dp.add_field(chr, 0);
	dp.add_field(pos, 1);
	dp.add_field(ref, 2);
	dp.add_field(alt, 3);
	
	dp.add_matrix(V_mat, false, 4, n_hsq);
	
	dp.parse_file(anchor_path, global_opts::global_region);
	dp.clear();
	
	if( chr.size() != n_snps ){
		std::cerr << "Number of variants in gvar file does not match genotype matrix dimensions." << "\n";
		abort();
	}
	
	V_mat = (V_mat * (getPredParams(hsq_vals).transpose())).eval();
	
	return;
}

void save_V_anchor_points( bcf_srs_t*& sr, bcf_hdr_t*& hdr, genotype_data& g_data, table& c_data, Eigen::SparseMatrix<double>& L, Eigen::VectorXd& GRM_lambda )
{
	
	Eigen::SparseMatrix<double> Q;
	Eigen::VectorXd Q_lambda;

	subset_eigen(L, GRM_lambda, Q, Q_lambda);

	Eigen::MatrixXd& C = c_data.data_matrix;
	
	std::cerr << "Calculating partial rotations ...\n";
	Eigen::MatrixXd QtY, CtY, QtC, QtG, CtC, CtC_i;

	QtG = (Q.transpose() * g_data.genotypes).eval();
	
	QtC = (Q.transpose() * C).eval();
	CtC = (C.transpose() * C).eval();
	CtC_i = (CtC.inverse()).eval();
	
	std::cerr << "Rotating expression and covariates ... ";
	
	
	std::cerr << "Done.\n";

	std::vector<double> hsq_vals{0.0, 0.5, 1.0};
	
	Eigen::MatrixXd V_mat, X;
	
	calculate_V_anchor_points(V_mat, g_data, C, hsq_vals, CtC, CtC_i, QtG, QtC, Q_lambda);
	
	return;
}
