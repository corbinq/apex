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

#include "fitModels.hpp"

static const double min_hsq_thresh = 1e-5;
// static const int n_fa_iter = 3;
static const double adj_uniq = 1e-3;


void calc_tSVD(const Eigen::MatrixXd& Y, Eigen::MatrixXd& U, Eigen::VectorXd& d, Eigen::MatrixXd& V, const int& nv){

	Eigen::MatrixXd R = Y * Y.transpose();

	Spectra::DenseSymMatProd<double> op(R);
	
	int n2 = 2* nv;
	if( n2 > R.cols() ){
		n2 = R.cols();
	}
	
	Spectra::SymEigsSolver< double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double> > eigs(&op, nv, n2);

	eigs.init();
	int nconv = eigs.compute();
	
	d = eigs.eigenvalues().cwiseAbs().cwiseSqrt();
	U = eigs.eigenvectors();
	V = Y.transpose() * U * d.asDiagonal().inverse();
	
	return;
}


void calc_left_tSVD(const Eigen::MatrixXd& Y, Eigen::MatrixXd& U, Eigen::VectorXd& d, const int& nv){

	Eigen::MatrixXd R = Y * Y.transpose();

	Spectra::DenseSymMatProd<double> op(R);
	
	int n2 = 2* nv;
	if( n2 > R.cols() ){
		n2 = R.cols();
	}
	
	Spectra::SymEigsSolver< double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double> > eigs(&op, nv, n2);

	eigs.init();
	int nconv = eigs.compute();
	
	d = eigs.eigenvalues().cwiseAbs().cwiseSqrt();
	U = eigs.eigenvectors();
	
	return;
}


void calc_eGRM_PCs(Eigen::MatrixXd& ePCs, Eigen::VectorXd& lam, const Eigen::MatrixXd& Y, const int& n_eigs){
	
	// Eigen::MatrixXd YT = scale_and_center(Y.transpose());
	double m = Y.cols();
	double n = Y.rows();
	
	double adj = n/((n-1)*(m - 1));
	
	if(  global_opts::n_fa_iter == 0  ){
	
		// Makes trace(R) = n if Y is scaled.
		Eigen::MatrixXd R = Y * Y.transpose() * adj;
		
		// Ensure that trace(R) = n.
		adj = R.cols() / ( R.trace() );
		R *= adj;
		
		Spectra::DenseSymMatProd<double> op(R);
		
		int n2 = 2* n_eigs;
		if( n2 > R.cols() ){
			n2 = R.cols();
		}
		
		Spectra::SymEigsSolver< double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double> > eigs(&op, n_eigs, n2);

		eigs.init();
		int nconv = eigs.compute();
		
		lam = eigs.eigenvalues();
		ePCs = eigs.eigenvectors();
		
		// std::cout << "Values:\n";
		// std::cout << lam << "\n";
		// std::cout << "Vectors:\n";
		// std::cout << ePCs << "\n";
		
	}else{
		Eigen::MatrixXd Y_ = Y;
		scale_and_center(Y_);
		// Y_ *= std::sqrt(adj);
		
		//Eigen::MatrixXd Y_adj = Y_;
		
		Eigen::VectorXd s2_;
		
		s2_ = Y_.cwiseAbs2().colwise().mean();
		// s2_.array() += adj_uniq;
		
		// std::cout << "RESIDUAL VARIANCES 0:\n";
		// std::cout << s2_.mean() << "\n";
		// std::cout << s2_.minCoeff() << "\n";
		// std::cout << s2_.maxCoeff() << "\n";
	
		Eigen::MatrixXd U;
		Eigen::VectorXd d;
		// Eigen::MatrixXd V;
		
		for (int ii = 0; ii < global_opts::n_fa_iter; ii++){
			
			Eigen::MatrixXd Y_scaled = Y_ * s2_.cwiseAbs().cwiseSqrt().asDiagonal().inverse();
			
			//calc_tSVD(Y_scaled, U, d, V, n_eigs);
			calc_left_tSVD(Y_scaled, U, d, n_eigs);

			//Y_adj = U.transpose()Y_;
			//Y_adj = (Y_ - (U * d.asDiagonal() * V.transpose() * s2_.cwiseSqrt().asDiagonal())).eval();
			
			// s2_ = Y_adj.cwiseAbs2().colwise().mean();
			// s2_.array() += adj_uniq;
			
			for( int jj = 0; jj < Y_.cols(); jj++){
				Eigen::VectorXd Y_hat = U.transpose() * Y_.col(jj);
				Y_hat = (U * Y_hat).eval();
				double sigma2 = ((Y_.col(jj) - Y_hat).squaredNorm()/n);
				
				s2_(jj) = global_opts::fa_p * global_opts::fa_tau + ( 1 - global_opts::fa_p ) * sigma2;
				
			}
				
			// std::cout << "RESIDUAL VARIANCES " << ii << ":\n";
			// std::cout << s2_.mean() << "\n";
			// std::cout << s2_.minCoeff() << "\n";
			// std::cout << s2_.maxCoeff() << "\n";
			
		}
		
		lam = d.cwiseAbs2();
		lam.array() *= adj;
		ePCs = U;
	}
	
	return;
}


DiagonalXd calc_Vi(const double& phi, const Eigen::VectorXd& lambdas){
	Eigen::VectorXd vals = 1.00 + phi*lambdas.array();
	return vals.asDiagonal().inverse();
}

DiagonalXd calc_Psi(const double& phi, const Eigen::VectorXd& lambdas){
	
	double hsq = phi/(1.00 + phi);
	
	if( hsq <= min_hsq_thresh ){
		return Eigen::VectorXd::Zero(lambdas.size()).asDiagonal();
	}
	Eigen::VectorXd numer = hsq * lambdas;
	Eigen::VectorXd denom = hsq * lambdas;
	denom.array() += 1.00;
	return (numer.cwiseQuotient(denom)).asDiagonal();	
}

DiagonalXd calc_Psi_low_rank(const double& phi, const Eigen::VectorXd& lambdas){
	
	// double hsq = phi/(1.00 + phi);
	
	// if( hsq <= min_hsq_thresh ){
		// return Eigen::VectorXd::Zero(lambdas.size()).asDiagonal();
	// }
	Eigen::VectorXd numer = phi * lambdas;
	Eigen::VectorXd denom = lambdas;
	// denom.array() -= 1.00;
	denom.array() *= phi;
	denom.array() += 1.00;
	return (numer.cwiseQuotient(denom)).asDiagonal();	
}

Eigen::MatrixXd getVBeta(const std::vector<double>& hsq_v, const std::vector<double>& phi_v, const int& p){

	int q = hsq_v.size();

	Eigen::MatrixXd Beta(p, q);
	
	for(int i = 0; i < p; i++){
		for( int j = 0; j < q; j++){
			if( i == 0 ){
				Beta(i, j) = 1.00/(1.00 + phi_v[j]);
			}else{
				Beta(i, j) = std::pow(hsq_v[j], (double) i)/(1.00 + phi_v[j]);
			}
		}
	}
	return Beta;
}

Eigen::MatrixXd getPredParams( const std::vector<double>& vals ){
	int n = vals.size();
	Eigen::MatrixXd X(n,n);
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			if ( i == 0 && j == 0 ){
				X(i,j) = 1.0;
			}else{
				X(i,j) = std::pow(vals[i], (double) j );
			}
		}
	}
	return (X.transpose() * X).ldlt().solve(X.transpose());
}

Eigen::VectorXd predV( const Eigen::MatrixXd& vv, const double& hsq ){
	Eigen::VectorXd out = vv.col(0);
	
	Eigen::VectorXd Beta(vv.cols());
	
	Beta(0) = 1.00;
	
	for(int i = 1; i < vv.cols(); i++){
		Beta(i) = std::pow(hsq, (double) i);
	}
	return vv * Beta;
}

void block_diag_eigendecomp( Eigen::SparseMatrix<double>& M, std::vector<int>& block_sizes, Eigen::SparseMatrix<double>& eig_vectors, Eigen::VectorXd& eig_values ){
	
	using td = Eigen::Triplet<double>;

	int n = M.cols();
	int n_blocks = block_sizes.size();
	
	if( M.rows() != n ){
		std::cerr << "FATAL ERROR: Cannot decompose non-square block matrix!\n";
		abort();
	}

	int nvals = n;
	for( const int& nb : block_sizes){
		nvals += ( nb*(nb - 1) ) ;
	}

	std::vector<td> triplets;
	triplets.reserve(nvals);
	
	eig_values.resize(n);
	eig_vectors.resize(n,n);
	
	int offset = 0;
	
	// Perform decompositions 
	
	std::cerr << "Decomposing block";
	
	for( int i = 0; i < n_blocks; i++ ){
		
		const int& n_i = block_sizes[i];
		
		std::cerr << "\rDecomposing block " << i << " (" << n_i << ") out of " << n_blocks << " ... ";
		
		Eigen::MatrixXd M_i = Eigen::MatrixXd(M.topLeftCorner(n_i + offset, n_i + offset).bottomRightCorner(n_i, n_i));
		// M_i.makeCompressed();

		// std::cerr << "GRM eigendecomposition ... ";
		
		Eigen::SelfAdjointEigenSolver <Eigen::MatrixXd> Eig_i(M_i);
		
		if (Eig_i.info() != Eigen::Success){
			std::cerr << "FATAL ERROR: block decomposition failed!\n";
			abort();
		}
		
		for( int ii = 0; ii < n_i; ii++ ){
			
			eig_values.coeffRef( ii + offset ) = Eig_i.eigenvalues().coeffRef( ii );
			
			for( int jj = 0; jj < n_i; jj++ ){
				
				triplets.push_back(td(ii + offset,jj + offset, Eig_i.eigenvectors().coeffRef(ii,jj) ));
			}
		}
		
		offset += n_i;
	}
	std::cerr << "Done\n.";
	
	// Eigenvalues for unrelated block (identity matrix)
	for(int ii = offset; ii < n; ii++){
		eig_values.coeffRef(ii) = M.coeffRef(ii,ii);
		triplets.push_back(td(ii, ii, 1.00 ));
	}
	
	// Clear GRM, as we are now done with it. 
	M.resize(0,0);
	
	// Construct combined matrix of eigenvectors

	std::cerr << "Creating eigenvector matrix ... ";
	eig_vectors.setFromTriplets(triplets.begin(), triplets.end());
	std::cerr << "Done.\n";
	
	return;
}


void GRM_decomp( Eigen::SparseMatrix<double>& GRM, const std::vector<std::vector<int>>& relateds, Eigen::SparseMatrix<double>& L, Eigen::VectorXd& L_lambda, Eigen::SparseMatrix<double>& Q, Eigen::VectorXd& Q_lambda ){
	
	
	double delta_thresh = 1e-5;
	double drop0_thresh = 2e-6;

	int n = GRM.rows();

	int max_block_size = 0;
	std::vector<int> block_sizes;

	std::vector<bool> is_unrelated(n, true);
	
	std::vector<int> unrelateds;
	PermutXd Tr(n); 
	
	// Identify unrelated individuals, and record block sizes. 
	for( const auto& related_block_i : relateds ){
		int size_i = related_block_i.size();
		if( size_i > max_block_size ){
			max_block_size = size_i;
		}
		block_sizes.push_back( size_i );
		
		for( const auto& ii : related_block_i ){
			is_unrelated[ii] = false;
		}
	}
	// Collect ordered indices of unrelated individuals. 
	for( int j = 0; j < n; j++ ){
		if( is_unrelated[j] ){
			unrelateds.push_back(j);
		}
	}
	is_unrelated.clear();
	
	int i = 0;
	
	for( const auto& related_block_i : relateds ){
		for( const int& ii : related_block_i ){
			Tr.indices()[ii] = i;
			i++;
		}
	}
	for( const int& ii : unrelateds ){
		Tr.indices()[ii] = i;
		i++;
	}

	int nind = unrelateds.size();
	int nrel = n - nind;
	int n_blocks = relateds.size();
	
	std::cerr << "Found " << n_blocks << " related blocks with " << nrel << " total individuals (largest block = " << max_block_size << ") ... \n";
	
	std::cerr << "Reordering related GRM blocks ... ";
	
	// A better way to do this would be  
	//    grm = grm.twistedBy(Tr);
	// Oddly, this fails on my system. 
	
	GRM = (Tr * GRM).eval();
	GRM = (GRM * Tr.transpose()).eval();
	
	GRM.makeCompressed();

	std::cerr << "Done.\n";

	// Eigen::SparseMatrix<double> GRM_rel = GRM.topLeftCorner(nrel,nrel).bottomRightCorner(nrel,nrel);
	// GRM_rel.makeCompressed();

	std::cerr << "GRM eigendecomposition ... ";
	
	block_diag_eigendecomp( GRM, block_sizes, L, L_lambda );
	
	std::cerr << "Done.\n";
	
	// Eigen::SelfAdjointEigenSolver <Eigen::SparseMatrix<double>> GRM_eig(GRM_rel);
	
	// if (GRM_eig.info() != Eigen::Success){
		// std::cerr << "FATAL ERROR: GRM decomposition failed!\n";
		// abort();
	// }
	// std::cerr << "Done.\n";
	
	// L_lambda = GRM_eig.eigenvalues();
	
	// L_lambda.conservativeResize(nrel + nind);
	// for(int i = nrel; i < nrel + nind; i++){
		// L_lambda(i) = GRM.coeffRef(i,i);
	// }
	
	
	// L = GRM_eig.eigenvectors().sparseView(drop0_thresh, 1.0 - std::numeric_limits<double>::epsilon());
	
	// L.conservativeResize(nrel + nind,nrel + nind);
	
	// for(int i = nrel; i < nrel + nind; i++){
		// L.coeffRef(i,i) = 1.00;
	// }
	
	L.makeCompressed();
	
	L = (Tr.transpose() * L).eval();
	
	L.makeCompressed();
	
	std::vector<int> kp;
	for( int i = 0; i < L_lambda.size(); i++ ){
		if( std::abs( L_lambda(i) - 1 ) > delta_thresh ){
			kp.push_back(i);
		}
	}
	
	// Permutation matrix to extract selected eigenvectors.
	using td = Eigen::Triplet<double>;
	Eigen::SparseMatrix<double> PC_sel(L_lambda.size(),kp.size());
	std::vector<td> PC_trips;
	
	// Selected eigenvalues.
	Q_lambda = Eigen::VectorXd(kp.size());
	
	i = 0;
	for( const int& k : kp){
		PC_trips.push_back(td(k,i,1.00));
		Q_lambda(i) =  L_lambda(k) - 1.00;
		i++;
	}
	PC_sel.setFromTriplets(PC_trips.begin(), PC_trips.end());
	
	std::cerr << "Selected " << kp.size() << " eigenvectors.  ";
	

	Q = (L * PC_sel).eval();
	Q.makeCompressed();
	
	std::cerr << "Done.\n";
	
	// GRM = (Tr.transpose() * GRM).eval();
	// GRM = (GRM * Tr).eval();
	
	// GRM.makeCompressed();

	return;
}


void GRM_decomp( Eigen::SparseMatrix<double>& GRM, const std::vector<std::vector<int>>& relateds, Eigen::SparseMatrix<double>& L, Eigen::VectorXd& L_lambda ){
	
	
	double delta_thresh = 1e-5;
	double drop0_thresh = 2e-6;

	int n = GRM.rows();

	int max_block_size = 0;
	std::vector<int> block_sizes;

	std::vector<bool> is_unrelated(n, true);
	
	std::vector<int> unrelateds;
	PermutXd Tr(n); 
	
	// Identify unrelated individuals, and record block sizes. 
	for( const auto& related_block_i : relateds ){
		int size_i = related_block_i.size();
		if( size_i > max_block_size ){
			max_block_size = size_i;
		}
		block_sizes.push_back( size_i );
		
		for( const auto& ii : related_block_i ){
			is_unrelated[ii] = false;
		}
	}
	// Collect ordered indices of unrelated individuals. 
	for( int j = 0; j < n; j++ ){
		if( is_unrelated[j] ){
			unrelateds.push_back(j);
		}
	}
	is_unrelated.clear();
	
	int i = 0;
	
	for( const auto& related_block_i : relateds ){
		for( const int& ii : related_block_i ){
			Tr.indices()[ii] = i;
			i++;
		}
	}
	for( const int& ii : unrelateds ){
		Tr.indices()[ii] = i;
		i++;
	}

	int nind = unrelateds.size();
	int nrel = n - nind;
	int n_blocks = relateds.size();
	
	std::cerr << "Found " << n_blocks << " related blocks with " << nrel << " total individuals (largest block = " << max_block_size << ") ... \n";
	
	std::cerr << "Reordering related GRM blocks ... ";
	
	// A better way to do this would be  
	//    grm = grm.twistedBy(Tr);
	// Oddly, this fails on my system. 
	
	GRM = (Tr * GRM).eval();
	GRM = (GRM * Tr.transpose()).eval();
	
	GRM.makeCompressed();

	std::cerr << "Done.\n";

	// Eigen::SparseMatrix<double> GRM_rel = GRM.topLeftCorner(nrel,nrel).bottomRightCorner(nrel,nrel);
	// GRM_rel.makeCompressed();

	std::cerr << "GRM eigendecomposition ... ";
	
	block_diag_eigendecomp( GRM, block_sizes, L, L_lambda );
	
	std::cerr << "reordering ... ";
	
	L.makeCompressed();
	
	L = (Tr.transpose() * L).eval();
	
	std::cerr << "Done.\n";
	
	return;
}


void subset_eigen( Eigen::SparseMatrix<double>& L, Eigen::VectorXd& L_lambda, Eigen::SparseMatrix<double>& Q, Eigen::VectorXd& Q_lambda ){
	
	double delta_thresh = 1e-5;

	std::vector<int> kp;
	for( int i = 0; i < L_lambda.size(); i++ ){
		if( std::abs( L_lambda(i) - 1 ) > delta_thresh ){
			kp.push_back(i);
		}
	}
	
	// Permutation matrix to extract selected eigenvectors.
	using td = Eigen::Triplet<double>;
	Eigen::SparseMatrix<double> PC_sel(L_lambda.size(),kp.size());
	std::vector<td> PC_trips;
	
	// Selected eigenvalues.
	Q_lambda = Eigen::VectorXd(kp.size());
	
	int i = 0;
	for( const int& k : kp){
		PC_trips.push_back(td(k,i,1.00));
		Q_lambda(i) =  L_lambda(k) - 1.00;
		i++;
	}
	PC_sel.setFromTriplets(PC_trips.begin(), PC_trips.end());
	
	std::cerr << "Selected " << kp.size() << " eigenvectors.  ";

	Q = (L * PC_sel).eval();
	Q.makeCompressed();
	
	return;
}


// void dense_GRM_decomp( Eigen::MatrixXd& GRM, const std::vector<std::vector<int>>& relateds, Eigen::MatrixXd& L, Eigen::VectorXd& L_lambda, Eigen::MatrixXd& Q, Eigen::VectorXd& Q_lambda ){
	
	// double delta_thresh = 1e-5;
	// double drop0_thresh = 2e-6;

	// std::vector<int> unrelateds;
	// PermutXd Tr(GRM.rows()); 
	
	// int i = 0;
	// for( int j = 0; j < GRM.cols(); j++ ){
		// if( i < relateds.size() ){
			// if( relateds[i] == j ){
				// i++;
			// }else{
				// unrelateds.push_back(j);
			// }
		// }else{
			// unrelateds.push_back(j);
		// }
	// }
	// i = 0;
	// for( const int& ii : relateds ){
		// Tr.indices()[ii] = i;
		// i++;
	// }
	// for( const int& ii : unrelateds ){
		// Tr.indices()[ii] = i;
		// i++;
	// }

	// int nrel = relateds.size();
	// int nind = unrelateds.size();
	
	// std::cerr << "Found "<< nrel << " related individuals ... \n";
	
	// std::cerr << "Reordering related GRM blocks ... ";
	
	// // A better way to do this would be  
	// //    grm = grm.twistedBy(Tr);
	// // Oddly, this fails on my system. 
	
	// GRM = (Tr * GRM).eval();
	// GRM = (GRM * Tr.transpose()).eval();
	
	// std::cerr << "Done.\n";

	// std::cerr << "GRM eigendecomposition ... ";
	
	// Eigen::SelfAdjointEigenSolver <Eigen::MatrixXd> GRM_eig(GRM.topLeftCorner(nrel,nrel));
	
	// if (GRM_eig.info() != Eigen::Success){
		// std::cerr << "FATAL ERROR: GRM decomposition failed!\n";
		// abort();
	// }
	// std::cerr << "Done.\n";
	
	// L_lambda = GRM_eig.eigenvalues();
	
	// L_lambda.conservativeResize(nrel + nind);
	// for(int i = nrel; i < nrel + nind; i++){
		// L_lambda(i) = GRM.coeffRef(i,i);
	// }
	
	
	// L = GRM_eig.eigenvectors();
	
	// L.conservativeResize(nrel + nind,nrel + nind);
	
	// for(int i = nrel; i < nrel + nind; i++){
		// L.coeffRef(i,i) = 1.00;
	// }
	
	// L = (Tr.transpose() * L).eval();
	
	// std::vector<int> kp;
	// for( int i = 0; i < L_lambda.size(); i++ ){
		// if( std::abs( L_lambda(i) - 1 ) > delta_thresh ){
			// kp.push_back(i);
		// }
	// }
	
	// // Permutation matrix to extract selected eigenvectors.
	// using td = Eigen::Triplet<double>;
	// Eigen::MatrixXd PC_sel(L_lambda.size(),kp.size());
	
	// // Selected eigenvalues.
	// Q_lambda = Eigen::VectorXd(kp.size());
	
	// i = 0;
	// for( const int& k : kp){
		// PC_sel(k,i) = 1.00;
		// Q_lambda(i) =  L_lambda(k) - 1.00;
		// i++;
	// }
	
	// std::cerr << "Selected " << kp.size() << " eigenvectors.\n";
	
	// Q = (L * PC_sel).eval();
	
	// return;
// }


void meta_svar_sumstat::condition_on_het(const int& k){
	kept_snps.push_back(k);
	for( int s = 0; s < ss.size(); s++ ){
		//Eigen::MatrixXd cov_mat = vg.Covar_perStudy(s, k);
		//cout << cov_mat << "\n";
		ss[s].condition_on( std::vector<int>(1, k), vg.Covar_perStudy(s, k) );
		// std::cout << "Conditioned " << s << "\n";
	}
	update_meta_ss();
}

lm_output lm_from_sumstats( const Eigen::VectorXd& U, const Eigen::VectorXd& V, const double& n, const double& df_0, const double& stdev, const Eigen::VectorXd& U_0, const Eigen::MatrixXd& J_0, const Eigen::MatrixXd& Cov, const bool& check_filter, const std::vector<bool>& exclude ){
	
	lm_output out;
	
	Eigen::MatrixXd Ji;
	Eigen::VectorXd JiU_0;
	Eigen::VectorXd C_JiU_0;
	
	double m_0 = U_0.size();
	
	double SSE_0 = df_0;
	
	double df = df_0 - m_0 - 1;
	
	int i_m = 0;
	
	if( m_0 > 0 ){
		
		Ji = J_0.inverse();

		JiU_0 = Ji * U_0;

		C_JiU_0 = Cov * JiU_0;

		SSE_0 = df_0 - U_0.dot(JiU_0);
	}
	
	for( int i = 0; i < U.size(); ++i){
		bool skip = false;
		if( exclude.size() > 0 ) skip = exclude[i];
		if( skip && check_filter ){
			// std::cerr << "\nskip\n\n";
			out.push_back(-99, -99, -99, NAN);
		}else{
			double SCORE, VARSC, SSE_i;
			SSE_i = SSE_0;
			if( m_0 > 0 ){
				SCORE = U(i) - C_JiU_0(i);
				VARSC = V(i) - Cov.row(i)*Ji*Cov.row(i).transpose();
			}else{
				SCORE = U(i);
				VARSC = V(i);
			}
			double VIF = calc_vif(V(i), VARSC);
			if( (VARSC > 0.0 && ((VARSC/V(i)) > (1.0 - global_opts::RSQ_PRUNE) && V(i) > 0.0)) || !check_filter ){
				double beta = stdev * SCORE/VARSC;
				
				double se = stdev * std::sqrt(SSE_i - SCORE*SCORE/VARSC) / std::sqrt(df*VARSC);
				long double pval = -99;
        double log_pval = NAN;

				double PVAR = SCORE*SCORE/VARSC;
				double STAT = df*PVAR/(SSE_i - PVAR);
				
				if( SSE_i - PVAR > 0 && STAT > 0 && PVAR > 0 ){
					log_pval = rmath::pf(STAT, 1.0, df, false, true); // calculate p-value with rmath pf in log scale
          pval = exp(log_pval);
				}

				out.push_back(beta, se, pval, log_pval, VIF);
			}else{
				// std::cerr << "\nWARNING: RSQ_VIF = " << (VARSC/V(i)) << ", VARSC = "<< VARSC << "\n\n";
				out.push_back(-99, -99, -99, NAN, VIF);
			}
		}
	}
	
	return out;
}

template <typename T>
int which_min( const std::vector<T>& p, bool gt0 ){
	T mp = p[0];
	int wm = -1;
	for(int i = 0; i < p.size(); ++i){
		if( gt0 && mp < 0 ){
			mp = p[i];
		}
		if( p[i] <= mp ){
			if( !(gt0 && p[i] <= 0) ){
				wm = i;
				mp = p[i];
			}
		}
	}
	return wm;
}

template <typename T>
int which_min(const std::vector<T>& v){
  T mp = std::numeric_limits<T>::max();
  int wm = -1;
  for(int i = 0; i < v.size(); ++i) {
    const T& val = v[i];
    if (std::isfinite(val)) {
      if (val < mp) {
        mp = val;
        wm = i;
      }
    }
  }
  return wm;
}

void forward_lm::check_joint_pvalues(int& index_of_largest_pvalue, double& largest_log_pvalue, const Eigen::VectorXd& U, const Eigen::VectorXd& V, const Eigen::VectorXd& U_0, const Eigen::MatrixXd& J_0, const Eigen::MatrixXd& Cov, const double& n, const double& m){
	
	index_of_largest_pvalue = -1;
  largest_log_pvalue = -std::numeric_limits<double>::infinity();
	
	int k_i = 0;
	std::vector<int> kept_snps = seq_int(keep.size());
	for(const int k : keep ){
		std::vector<int> kept_snps_not_k = kept_snps;
		kept_snps_not_k.erase( kept_snps_not_k.begin() + k_i );
	
		double current_log_pvalue;

		Eigen::VectorXd U_k = U(std::vector<int>(1,k));
		Eigen::VectorXd V_k = V(std::vector<int>(1,k));
			
		Eigen::VectorXd U_0k = U_0(kept_snps_not_k); 
		Eigen::MatrixXd J_0k = J_0(kept_snps_not_k, kept_snps_not_k);
    Eigen::MatrixXd Cov_k = Cov(std::vector<int>(1,k), kept_snps_not_k);

		lm_output reg_k = lm_from_sumstats(U_k, V_k, n, m, 1.00, U_0k, J_0k, Cov_k, false);

    current_log_pvalue = reg_k.log_pval[0];

		if(current_log_pvalue > largest_log_pvalue || isnan(current_log_pvalue)){
      largest_log_pvalue = current_log_pvalue;
			index_of_largest_pvalue = k_i;
		}
		
		k_i++;
	}
}

// Used by conditional_analysis()
forward_lm::forward_lm(const Eigen::VectorXd& U, const Eigen::VectorXd& V, const double& n, const double& m, const double& stdev, vcov_getter& vget, double pval_thresh, const std::vector<double>& weights )
{

	Eigen::VectorXd U_0 = Eigen::VectorXd(0);
	Eigen::MatrixXd J_0 = Eigen::MatrixXd(0,0);
	Eigen::MatrixXd Cov = Eigen::VectorXd(0);
	
	if( U.size() != V.size() ){
		std::cerr << "U.size() != V.size()" << "\n";
		exit(1);
	}
	
	if( U.size() <= 0 ){
		std::cerr << "U.size() <= 0" << "\n";
		exit(1);
	}
	
	int n_var = U.size();
	int nk = 0;
	
	std::vector<bool> excl(n_var,false);
	
	lm_output reg0;
	
	double alpha_thresh = global_opts::LM_ALPHA;
  double log_alpha = log(alpha_thresh);
	
	int total_steps = 0;
	
	while( 1 )
	{
		
		int steps_taken = 0;
		
		lm_output reg = lm_from_sumstats(U, V, n, m, stdev, U_0, J_0, Cov, true, excl);
		// std::cerr << "Fit model.\n";

		if( nk == 0 && reg0.beta.size() == 0 ){
			reg0 = reg;
		}

		int wk = which_min(reg.log_pval);
		double adj_pval = ACAT_non_missing(reg.pval, weights);
    double log_adj_pval = log(adj_pval);

		//cout << nk << ":" << wk << "\t" << reg0.beta[wk]  << "\t" << reg0.se[wk]  << "\t" << reg.pval[wk]  << "\t" << adj_pval << "\n";
		
		double log_pval_check = NAN;
		
		if( wk >= 0 ){
			if( global_opts::step_marginal ){
        log_pval_check = reg.log_pval[wk];
			}else{
        log_pval_check = log_adj_pval;
			}
		}
		
		// -----------------------------------
		// Forward step. 
		// -----------------------------------
		if(  (!isnan(log_pval_check) && log_pval_check < log_alpha) || (keep.size() == 0 && wk >= 0 )  )
		{
			keep.push_back(wk);
			// excl[wk] = true;
			
			beta_0.push_back(reg0.beta[wk]);
			se_0.push_back(reg0.se[wk]);
			pval_0.push_back(reg0.pval[wk]);
			pval_seq.push_back(reg.pval[wk]);
			pval_adj.push_back(adj_pval);
      log_pval_0.push_back(reg0.log_pval[wk]);
      log_pval_seq.push_back(reg.log_pval[wk]);
      log_pval_adj.push_back(log_adj_pval);
			
			nk++;
			
			if( nk == 1 ){
				U_0 = Eigen::VectorXd(1);
				Cov = Eigen::MatrixXd(n_var,1);

			}else{
				U_0.conservativeResize(nk);
				Cov.conservativeResize(Eigen::NoChange, nk);
			}
			

			//get_var(keep, J_0);
			J_0 = vget.Var(keep);

			U_0(nk-1) = U(wk);
			
			if( log_pval_check <= log_alpha || keep.size() < 2 ){
				// Eigen::VectorXd new_cov; 
				// get_cov(wk, new_cov);
				Eigen::VectorXd new_cov = vget.Covar(wk);


				if( new_cov.size() != Cov.rows() )
				{
					std::cerr << "new_cov.size() != Cov.rows()" << "\n";
					exit(1);
				}
				
				std::vector<int> b_list;
				std::vector<double> r_list;
				
				for(int i = 0; i < n_var; ++i)
				{
					
					Cov(i, nk-1) = new_cov(i);
					
					if( global_opts::RSQ_BUDDY < 1.00 ){
						double corr = new_cov(i)/std::sqrt(V(i)*V(wk));
						
						// if( corr*corr > global_opts::RSQ_PRUNE ) excl[i] = true;
						
						if( corr* corr > global_opts::RSQ_BUDDY ){
							b_list.push_back(i);
							r_list.push_back(corr);
						}
					}
				}
				
				if( global_opts::RSQ_BUDDY < 1.00 ){
					buddy_list.push_back(b_list);
					corr_list.push_back(r_list);
				}
			}
			
			steps_taken++;
			
			step_history.push_back(step_record{total_steps, steps_taken, 1, wk, keep.size() });
			
			//cout << "forward: " << nk << " SNPs. Added SNP " << wk << "\n";
		}
		
		// -----------------------------------
		// Backward step. 
		// -----------------------------------
		if( global_opts::backward_thresh < 1.00 && nk > 1 ){
      double log_backthresh = log(global_opts::backward_thresh);
			double max_log_joint_pvalue = -std::numeric_limits<double>::infinity();
			int k_rm = -1;
			
			check_joint_pvalues(k_rm, max_log_joint_pvalue, U, V, U_0, J_0, Cov, n, m);
			//cout << max_joint_pvalue << "\n";
			
			if( (max_log_joint_pvalue > log_backthresh || isnan(max_log_joint_pvalue) ) && k_rm < nk - 1 && k_rm >= 0 ){
				
				std::vector<int> kept_snps_not_k_rm = seq_int(nk);
				kept_snps_not_k_rm.erase( kept_snps_not_k_rm.begin() + k_rm );
		
				U_0 = (U_0(kept_snps_not_k_rm)).eval(); 
				J_0 = (J_0(kept_snps_not_k_rm, kept_snps_not_k_rm)).eval(); 
				Cov = (Cov(seq_int(n_var), kept_snps_not_k_rm)).eval();
					
				beta_0.erase( beta_0.begin() + k_rm );
				se_0.erase( se_0.begin() + k_rm );

				pval_0.erase( pval_0.begin() + k_rm );
				pval_seq.erase( pval_seq.begin() + k_rm );
				pval_adj.erase( pval_adj.begin() + k_rm );

        log_pval_0.erase( log_pval_0.begin() + k_rm );
        log_pval_seq.erase( log_pval_seq.begin() + k_rm );
        log_pval_adj.erase( log_pval_adj.begin() + k_rm );

				
				if( global_opts::RSQ_BUDDY < 1.00 ){
					buddy_list.erase(buddy_list.begin() + k_rm);
					corr_list.erase(corr_list.begin() + k_rm);
				}
				
				// excl[keep[k_rm]] = false;
				
				int snp_dropped = keep[k_rm];
				
				keep.erase( keep.begin() + k_rm );
				
				nk--;
				
				// std::cout << "backward\n";
				steps_taken++;
				
				step_history.push_back(step_record{total_steps, steps_taken, -1, snp_dropped, keep.size() });
			}
		}
		// std::cout << steps_taken << "\n";
		
		if( steps_taken == 0 ){
			step_history.push_back(step_record{total_steps, steps_taken, 0, 0, keep.size() });
			
			break;
		}else if( keep.size() >= global_opts::max_signals ){
			step_history.push_back(step_record{total_steps, steps_taken, 0, 1, keep.size() });
			
			break;
		}
		total_steps++;
		if( total_steps > global_opts::max_steps ){
			step_history.push_back(step_record{total_steps, steps_taken, 0, 2, keep.size() });
			std::cerr << "\n\nWARNING: Exceeded max steps. Convergence failed.\n\n";
			break;
		}
	}

	int k_i = 0;
	std::vector<int> kept_snps = seq_int(nk);
	
	for( const int k : keep ){
		
		std::vector<int> kept_snps_not_k = kept_snps;
		kept_snps_not_k.erase( kept_snps_not_k.begin() + k_i );
	
		if( kept_snps_not_k.size() == 0 ){
			
			beta.push_back(reg0.beta[k]);
			se.push_back(reg0.se[k]);
			pval_joint.push_back(reg0.pval[k]);
      log_pval_joint.push_back(reg0.log_pval[k]);

		}else{

      /**
       * U_k is vector of length 1, containing score stat of variant k
       * V_k ^                    , containing variance of score stat of variant k
       * U_0k is vector of length nk, score statistics for all variants except k
       * J_0k is matrix of (nk-1) x (nk-1), covariance for all variants in model except k
       * Cov_k is matrix of 1 x (nk-1), covariance between variant k and all other variants
       */
			Eigen::VectorXd U_k = U(std::vector<int>(1,k));
			Eigen::VectorXd V_k = V(std::vector<int>(1,k));
				
			Eigen::VectorXd U_0k = U_0(kept_snps_not_k); 
			Eigen::MatrixXd J_0k = J_0(kept_snps_not_k, kept_snps_not_k); 
			Eigen::MatrixXd Cov_k = Cov(std::vector<int>(1,k), kept_snps_not_k);
			
			lm_output reg_k = lm_from_sumstats(U_k, V_k, n, m, stdev, U_0k, J_0k, Cov_k, false);
			
			beta.push_back(reg_k.beta[0]);
			se.push_back(reg_k.se[0]);
			pval_joint.push_back(reg_k.pval[0]);
      log_pval_joint.push_back(reg_k.log_pval[0]);
		}

		conditioned.push_back(kept_snps);
		k_i++;
		
	}
	
}

// to do: make parent generic class vcov_getter, where indiv_vcov_getter and sumstat_vcov_getter inherit from vcov_getter. This is an exact duplicate function. 

forward_lm::forward_lm(const Eigen::VectorXd& U, const Eigen::VectorXd& V, const double& n, const double& m, const double& stdev, indiv_vcov_getter& vget, double pval_thresh, const std::vector<double>& weights )
{

	Eigen::VectorXd U_0 = Eigen::VectorXd(0);
	Eigen::MatrixXd J_0 = Eigen::MatrixXd(0,0);
	Eigen::MatrixXd Cov = Eigen::VectorXd(0);
	
	if( U.size() != V.size() ){
		std::cerr << "U.size() != V.size()" << "\n";
		exit(1);
	}
	
	if( U.size() <= 0 ){
		std::cerr << "U.size() <= 0" << "\n";
		exit(1);
	}
	
	int n_var = U.size();
	int nk = 0;
	
	std::vector<bool> excl(n_var,false);
	
	lm_output reg0;
	
	double alpha_thresh = global_opts::LM_ALPHA;
	
	int total_steps = 0;
	
	while( 1 )
	{
		
		int steps_taken = 0;
		
		lm_output reg = lm_from_sumstats(U, V, n, m, stdev, U_0, J_0, Cov, true, excl);
		// std::cerr << "Fit model.\n";

		if( nk == 0 && reg0.beta.size() == 0 ){
			reg0 = reg;
		}

		int wk = which_min(reg.pval, true);
		double adj_pval = ACAT_non_missing(reg.pval, weights);
		
		//cout << nk << ":" << wk << "\t" << reg0.beta[wk]  << "\t" << reg0.se[wk]  << "\t" << reg.pval[wk]  << "\t" << adj_pval << "\n";
		
		long double pval_check = -99;
		
		if( wk >= 0 ){
			if( global_opts::step_marginal ){
				pval_check = reg.pval[wk];
			}else{
				pval_check = adj_pval;
			}
		}
		
		// -----------------------------------
		// Forward step. 
		// -----------------------------------
		if(  (pval_check >= 0 && pval_check < alpha_thresh) || (keep.size() == 0 && wk >= 0 )  )
		{
			keep.push_back(wk);
			// excl[wk] = true;
			
			beta_0.push_back(reg0.beta[wk]);
			se_0.push_back(reg0.se[wk]);
			pval_0.push_back(reg0.pval[wk]);
			pval_seq.push_back(reg.pval[wk]);
			pval_adj.push_back(adj_pval);
			
			nk++;
			
			if( nk == 1 ){
				U_0 = Eigen::VectorXd(1);
				Cov = Eigen::MatrixXd(n_var,1);

			}else{
				U_0.conservativeResize(nk);
				Cov.conservativeResize(Eigen::NoChange, nk);
			}
			

			//get_var(keep, J_0);
			J_0 = vget.Var(keep);

			U_0(nk-1) = U(wk);
			
			if( pval_check <= alpha_thresh || keep.size() < 2 ){
				// Eigen::VectorXd new_cov; 
				// get_cov(wk, new_cov);
				Eigen::VectorXd new_cov = vget.Covar(wk);


				if( new_cov.size() != Cov.rows() )
				{
					std::cerr << "new_cov.size() != Cov.rows()" << "\n";
					exit(1);
				}
				
				std::vector<int> b_list;
				std::vector<double> r_list;
				
				for(int i = 0; i < n_var; ++i)
				{
					
					Cov(i, nk-1) = new_cov(i);
					
					if( global_opts::RSQ_BUDDY < 1.00 ){
						double corr = new_cov(i)/std::sqrt(V(i)*V(wk));
						
						// if( corr*corr > global_opts::RSQ_PRUNE ) excl[i] = true;
						
						if( corr* corr > global_opts::RSQ_BUDDY ){
							b_list.push_back(i);
							r_list.push_back(corr);
						}
					}
				}
				
				if( global_opts::RSQ_BUDDY < 1.00 ){
					buddy_list.push_back(b_list);
					corr_list.push_back(r_list);
				}
			}
			
			steps_taken++;
			//cout << "forward: " << nk << " SNPs. Added SNP " << wk << "\n";
		}
		
		// -----------------------------------
		// Backward step. 
		// -----------------------------------
		if( global_opts::backward_thresh < 1.00 && nk > 1 ){

			double max_joint_pvalue = 0; 
			int k_rm = -1;
			
			check_joint_pvalues(k_rm, max_joint_pvalue, U, V, U_0, J_0, Cov, n, m);
			//cout << max_joint_pvalue << "\n";
			
			if( (max_joint_pvalue > global_opts::backward_thresh || max_joint_pvalue < 0 ) && k_rm < nk - 1 && k_rm >= 0 ){
				
				std::vector<int> kept_snps_not_k_rm = seq_int(nk);
				kept_snps_not_k_rm.erase( kept_snps_not_k_rm.begin() + k_rm );
		
				U_0 = (U_0(kept_snps_not_k_rm)).eval(); 
				J_0 = (J_0(kept_snps_not_k_rm, kept_snps_not_k_rm)).eval(); 
				Cov = (Cov(seq_int(n_var), kept_snps_not_k_rm)).eval();
					
				beta_0.erase( beta_0.begin() + k_rm );
				se_0.erase( se_0.begin() + k_rm );
				pval_0.erase( pval_0.begin() + k_rm );
				
				pval_seq.erase( pval_seq.begin() + k_rm );
				pval_adj.erase( pval_adj.begin() + k_rm );
				
				if( global_opts::RSQ_BUDDY < 1.00 ){
					buddy_list.erase(buddy_list.begin() + k_rm);
					corr_list.erase(corr_list.begin() + k_rm);
				}
				
				// excl[keep[k_rm]] = false;
				
				keep.erase( keep.begin() + k_rm );
				
				nk--;
				
				// std::cout << "backward\n";
				steps_taken++;
			}
		}
		// std::cout << steps_taken << "\n";
		
		if( steps_taken == 0 ){
			break;
		}else if( keep.size() >= global_opts::max_signals ){
			break;
		}
		total_steps++;
		if( total_steps > global_opts::max_steps ){
			std::cerr << "\n\nWARNING: Exceeded max steps. Convergence failed.\n\n";
			break;
		}
	}

	int k_i = 0;
	std::vector<int> kept_snps = seq_int(nk);
	
	for( const int k : keep ){
		
		std::vector<int> kept_snps_not_k = kept_snps;
		kept_snps_not_k.erase( kept_snps_not_k.begin() + k_i );
	
		if( kept_snps_not_k.size() == 0 ){
			
			beta.push_back(reg0.beta[k]);
			se.push_back(reg0.se[k]);
			pval_joint.push_back(reg0.pval[k]);

		}else{
			
			Eigen::VectorXd U_k = U(std::vector<int>(1,k));
			Eigen::VectorXd V_k = V(std::vector<int>(1,k));
				
			Eigen::VectorXd U_0k = U_0(kept_snps_not_k); 
			Eigen::MatrixXd J_0k = J_0(kept_snps_not_k, kept_snps_not_k);
      Eigen::MatrixXd Cov_k = Cov(std::vector<int>(1,k), kept_snps_not_k);

			lm_output reg_k = lm_from_sumstats(U_k, V_k, n, m, stdev, U_0k, J_0k, Cov_k, false);
			
			beta.push_back(reg_k.beta[0]);
			se.push_back(reg_k.se[0]);
			pval_joint.push_back(reg_k.pval[0]);
		}

		conditioned.push_back(kept_snps);
		k_i++;
		
	}
	
}

