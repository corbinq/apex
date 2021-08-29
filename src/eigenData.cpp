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

#include "eigenData.hpp"


void read_eigen(const std::string& file_name, Eigen::SparseMatrix<double>& eig_vec, Eigen::VectorXd& eig_val, const std::vector<std::string>& kp_ids){
	
	int n_samples, n_eigs, nnzero;
	
	// std::string file_name = file_prefix + ".eigen.gz";
	// file_name = std::regex_replace(file_name, std::regex(".eigen.gz.eigen.gz"), std::string(".eigen.gz"));
	
	// ----------------------------------------------
	// ----------------------------------------------
	// Read file 
	
	htsFile *htsf = hts_open(file_name.c_str(), "r");
	kstring_t str = {0,0,0};
	int n_fields;
	int *offsets;
	
	
	// -------------------------
	// 1) First line contains number of samples, number of eigenvectors, number non-zero. 
	
	hts_getline(htsf, KS_SEP_LINE, &str);
	if( !str.l ){
		std::cerr << "ERROR: eigen file is incomplete.\n";
		abort();
	}
	offsets = ksplit(&str, 0, &n_fields);
	
	// Number of samples
	n_samples = std::atoi((str.s + offsets[0]));
	// NUmber of eigenvectors (retained based on eigenvalues)
	n_eigs = std::atoi((str.s + offsets[1]));
	// Note: Get this with sp.nonZeros();
	nnzero = std::atoi((str.s + offsets[2]));
	
	if( n_samples != kp_ids.size() ){
		std::cerr << "ERROR: sample IDs in eigen file do not match in input.\n";
		abort();
	}
	
	// Resize the vector and matrix objects
	eig_vec.resize(n_samples, n_eigs);
	eig_val.resize(n_eigs);
	
	
	// -------------------------
	// 2) Second line contains number of samples and number of eigenvectors. 
	
	hts_getline(htsf, KS_SEP_LINE, &str);
	if( !str.l ){
		std::cerr << "ERROR: eigen file is incomplete.\n";
		abort();
	}
	offsets = ksplit(&str, 0, &n_fields);
	
	if( n_fields != kp_ids.size() ){
		std::cerr << "ERROR: sample IDs in eigen file do not match input.\n";
		abort();
	}
	
	for( int i = 0; i < n_fields; i++)
	{
		if(  std::string(str.s + offsets[i]) != kp_ids[i] ){
			std::cerr << "ERROR: sample IDs in eigen file do not match input.\n";
			abort();
		}
	}
	
	// -------------------------
	// 3) Third line contains eigenvalues. 
	
	
	hts_getline(htsf, KS_SEP_LINE, &str);
	if( !str.l ){
		std::cerr << "ERROR: eigen file is incomplete.\n";
		abort();
	}
	offsets = ksplit(&str, 0, &n_fields);
	
	if( n_fields != n_eigs ){
		std::cerr << "ERROR: Incorrect number of eigenvalues.\n";
		abort();
	}
	
	for( int i = 0; i < n_fields; i++)
	{
		eig_val.coeffRef(i) = std::atof(str.s + offsets[i]);
	}
	
	
	// -------------------------
	// 4) Remaining lines contain eigenvector triplets. 
	
	using td = Eigen::Triplet<double>;

	std::vector<td> triplets;
	triplets.reserve(nnzero);
	
	while( hts_getline(htsf, KS_SEP_LINE, &str) ){
		if( !str.l ) break;
		
		offsets = ksplit(&str, 0, &n_fields);
		
		if( n_fields != 3 ){
			std::cerr << "ERROR: Malformed eigen file.\n";
			// std::cerr << "Line " << triplets.size() << "\n";
			// std::cerr << "No. fields " << n_fields << "\n";
			abort();
		}
		
		int ii = std::atoi(str.s + offsets[0]);
		int jj = std::atoi(str.s + offsets[1]);
		double val = std::atof(str.s + offsets[2]);
		
		triplets.push_back(td(ii,jj,val));
	}
	
	eig_vec.setFromTriplets(triplets.begin(), triplets.end());
	
	// -------------------------
	// 5) Close input file. 
	
	ks_free(&str);
	hts_close(htsf);
}


void write_eigen(const std::string& file_prefix, Eigen::SparseMatrix<double>& eig_vec, Eigen::VectorXd& eig_val, const std::vector<std::string>& kp_ids){
	
	// std::cerr << "Start write eigen.\n"; // NOTE: DEBUG LINE
	
	// std::cerr << "Prefix: " << file_prefix << "\n";  // NOTE: DEBUG LINE
	
	std::string file_name = file_prefix + ".eigen.gz";
	
	
	// std::cerr << "Prefix: " << file_prefix << "\n";  // NOTE: DEBUG LINE
	
	// file_name = std::regex_replace(file_name, std::regex(".eigen.gz.eigen.gz"), std::string(".eigen.gz"));
	
	// std::cerr << "Processed file name.\n"; // NOTE: DEBUG LINE
	
	// std::cerr << file_name << "\n";  // NOTE: DEBUG LINE
	
	// ----------------------------------------------
	// ----------------------------------------------
	// Open output file  
	
	BGZF* out_file = bgzf_open(file_name.c_str(), "w");
	
	// std::cerr << "Opened file.\n"; // NOTE: DEBUG LINE
	
	// -------------------------
	// 1) First line contains number of samples, number of eigenvectors, number non-zero. 
	
	write_to_bgzf(std::to_string( kp_ids.size() ), out_file );
	write_to_bgzf("\t", out_file );
	write_to_bgzf(std::to_string( eig_val.size() ), out_file  );
	write_to_bgzf("\t", out_file );
	write_to_bgzf(std::to_string( eig_vec.nonZeros() ), out_file  );
	write_to_bgzf("\n", out_file );
	
	// std::cerr << "1\n"; // NOTE: DEBUG LINE
	
	// -------------------------
	// 2) Second line contains list of sample IDs
	
	write_to_bgzf(kp_ids[0], out_file);
	for( int i = 1; i < kp_ids.size(); i++ ){
		write_to_bgzf("\t", out_file );
		write_to_bgzf(kp_ids[i], out_file );
	}
	write_to_bgzf("\n", out_file );
	
	// std::cerr << "2\n"; // NOTE: DEBUG LINE
	
	// -------------------------
	// 3) Third line contains eigenvalues. 
	
	write_to_bgzf(std::to_string(eig_val.coeffRef(0)), out_file  );
	for( int i = 1; i < eig_val.size(); i++ ){
		write_to_bgzf("\t", out_file );
		write_to_bgzf(std::to_string(eig_val.coeffRef(i)), out_file  );
	}
	write_to_bgzf("\n", out_file);
	
	// std::cerr << "3\n"; // NOTE: DEBUG LINE
	
	// -------------------------
	// 4) Remaining lines contain eigenvector triplets. 
	
	for (int i = 0; i < eig_vec.outerSize(); i++)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator trip(eig_vec,i); trip; ++trip)
        {
			write_to_bgzf(std::to_string( trip.row() ), out_file  );
			write_to_bgzf("\t", out_file );
			write_to_bgzf(std::to_string( trip.col() ), out_file  );
			write_to_bgzf("\t", out_file );
			write_to_bgzf(std::to_string( trip.value() ), out_file  );
			write_to_bgzf("\n", out_file );
        }
    }
	
	// -------------------------
	// 5) Close the output file. 
	
	bgzf_close(out_file);
	
}



void update_blocks( int i, int j, std::vector<int>& cluster_ids, std::vector<std::set<int>>& clusters){
	if( i == j ){
		// A diagonal element, potentially block of size 1.
		// This does not affect block structure. 
		return;
	}
	
	int c_i = cluster_ids[i];
	int c_j = cluster_ids[j];
	
	if( c_i > c_j ){
		std::swap(c_i, c_j);
		std::swap(i, j);
	}
	
	if( c_i < 0 && c_j < 0 ){
		
		// Neither i nor j has been assigned to a block. 
		// Create a new block that comprises the pair. 
		
		cluster_ids[i] = clusters.size();
		cluster_ids[j] = clusters.size();
		clusters.push_back({i,j});
		
	}else if( c_i < 0 ){
		
		// Since c_i < c_j, only i might never have been assigned. 
		// If i was never assigned, then add it to j's block. 
		
		cluster_ids[i] = c_j;
		clusters[c_j].insert(i);
		
	}else if(c_i != c_j){
		
		// i and j were previously assigned to different blocks. 
		// We must merge these two blocks together. 
		
		// Add all members of j's block into i's block. 
		clusters[c_i].insert(clusters[c_j].begin(), clusters[c_j].end());
		
		// Set block ids of members of block c_j to c_i.
		for( const int& k : clusters[c_j] ){
			cluster_ids[k] = c_i;
		}
		
		// Delete the block c_j, as it has been merged into c_i.
		clusters.erase (clusters.begin() + c_j );
		
		// There is now one less block. 
		// For all blocks > c_j, decrement block id by 1. 
		for( int& k : cluster_ids ){
			if( k > c_j ){
				k -= 1;
			}
		}
		
	}
	return;
}


void read_sparse_GRM(const std::string& filename, Eigen::SparseMatrix<double>& GRM, const std::vector<std::string>& kp_ids, const double& r_scale, const int& r_col, std::vector<std::vector<int>>& related)
{
	int n = kp_ids.size();
	
	GRM = Eigen::SparseMatrix<double>(n,n);
	
	if( filename == "" ){
		GRM.setIdentity();
		return;
	}
	
	std::vector<int> cluster_ids;
	std::vector<std::set<int>> clusters;
	
	std::unordered_map<std::string, int> id_map;
	for( int i = 0; i < n; i++ ){
		id_map[ kp_ids[i] ] = i;
		cluster_ids.push_back(-1);
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
			
			update_blocks( ii, jj, cluster_ids, clusters);
			
			triplets.push_back(td(ii,jj,r_scale*val[i]));
			if( ii != jj ){
				triplets.push_back(td(jj,ii,r_scale*val[i]));
			}
		}
		
	}
	if( add_diag ){
		for(int i = 0; i < n; ++i) triplets.push_back(td(i,i,1.0));
	}
	
	related.clear();
	for( const auto& s : clusters ){
		related.push_back(std::vector<int>(s.begin(), s.end()));
	}
	
	GRM.setFromTriplets(triplets.begin(), triplets.end());
}
