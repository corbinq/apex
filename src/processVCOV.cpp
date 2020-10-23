#include "processVCOV.hpp"


double unflip(const int& f1, const int& f2){
	if( f1 != f2 ){
		return -1.0;
	}else{
		return 1.0;
	}
}

unsigned int pack_dp(const double& x, const double& m0, const double& m1){
	if( m0 > m1 ){
		return pack_dp(x, m1, m0);
	}else{
		if( 2.0 * m0 < VBIN_T_MAX ){
			return (unsigned int) round(x);
		}else{
			return (unsigned int) round(x/ceil((2.0*m0)/VBIN_T_MAX));
		}
	}
}

double unpack_dp(const double& x, const double& m0, const double& m1){
	if( m0 > m1 ){
		return unpack_dp(x, m1, m0);
	}else{
		if( 2.0 * m0 < VBIN_T_MAX ){
			return x;
		}else{
			return x*ceil((2.0*m0)/VBIN_T_MAX);
		}
	}
}

/*
vcov_bin_gz_out::vcov_bin_gz_out(std::string fn){
	open(fn);
}

void vcov_bin_gz_out::open(std::string fn){
	file_name = fn;
	bytes_c = 0;
	
	if( xz_mode ){
		buf_size = XZ_BUFFER_SIZE;
		std::string file_name_xz = file_name + ".xz";
		
		fp_xz = shrinkwrap::xz::ostream(file_name_xz.c_str());
	}else{
		buf_size = BUFFER_SIZE;
		std::string file_name_gz = file_name + ".gz";
		fp = bgzf_open(file_name_gz.c_str(), "w\0");

		bgzf_index_build_init(fp);
	}
}
*/

void vcov_bin_gz_out::add_to_queue(vbin_t val){
	if( queue.size() * sizeof(vbin_t) >= buf_size ){
		write_queue();
	}
	queue.push_back(val);
}

void vcov_bin_gz_out::write_queue(){

	size_t sz = sizeof(vbin_t);
	size_t nc = queue.size() * sz;
	
	if( xz_mode ){
		
		fp_xz.write(reinterpret_cast<const char*>(queue.data()), nc);
		fp_xz.flush();
		
	}else{
		
		void *buffer = malloc(buf_size);
		
		for( int ii = 0; ii < queue.size(); ii++ ){
			((vbin_t*)buffer)[ii] = queue[ii];
		}
		
		if (bgzf_write(fp, buffer, nc) < 0){
			std::cerr << "Fatal error: Couldn't write to " << file_name << ".gz\n";
			exit(1);
		} 
		
		free(buffer);
		
	}
	
	queue.clear();
}

void vcov_bin_gz_out::close(){
	
	write_queue();

	if( !xz_mode ){
		if ( bgzf_index_dump(fp, file_name.c_str(), ".gz.gzi") < 0 ){
			std::cerr << "Fatal error: Couldn't create " << file_name << ".gz.gzi\n";
			exit(1);
		}

		if (bgzf_close(fp) < 0){
			std::cerr << "Fatal error: Couldn't close " << file_name << ".gz\n";
			exit(1);
		}
	}
}

vcov_bin_file::vcov_bin_file(std::string fn){
	open(fn);
}

void vcov_bin_file::open(std::string file_name){
	
	is_xz = (file_name.substr( file_name.length() - 2 ) != "gz");
	
	if(is_xz){
		
		fp_xz.open(file_name);
		
	}else{
	
		fp = bgzf_open(file_name.c_str(), "r");
		
		if ( bgzf_index_load(fp, file_name.c_str(), ".gzi") < 0 ){
			std::cerr << "Error: Couldn't load index file for " << file_name << ".gz\n";
			exit(1);
		}
	}
}

void vcov_bin_file::close(){
	if(!is_xz){
		if (bgzf_close(fp) < 0){
			std::cerr << "Error: Couldn't close " << file_name << ".gz\n";
			exit(1);
		}
	}
}

std::vector<vbin_t> vcov_bin_file::get_data(long st, long n_vals){
	
	int c = 1;
    long start, end, size;
	
    start = st; 
	size = -1; end = -1;

	vbin_t value;
	int nval_buffer;
    
	//if (size >= 0) end = start + size;
	end = start + n_vals*sizeof(vbin_t);
	
	void* buffer = malloc(BUFFER_SIZE);

	std::vector<vbin_t> out_vec;

	if( !is_xz ){
		if ( bgzf_useek(fp, start, SEEK_SET) < 0 ){
			std::cerr << "Fatal error: Couldn't seek to byte " << start << " in " << file_name << "\n";
			exit(1);
		} 
	}

	while( start <= end && out_vec.size() < n_vals){
		
		int nbytes = BUFFER_SIZE;
		
		if (end >= 0){
			nbytes = (end - start > BUFFER_SIZE)? BUFFER_SIZE:(end - start);
		}
		
		if( !is_xz ){
			c = bgzf_read(fp, buffer, nbytes);
		}else{
			c = fp_xz.read(buffer, start, nbytes);
		}
		
		if (c == 0) break;
		if (c < 0){
			std::cerr << "Fatal error reading from VCOV file.\n";
			exit(1);
		}
		start += c;
		
		nval_buffer = (int)c/sizeof(vbin_t);
		
		for( int i = 0; i < nval_buffer; ++i ){
			if( out_vec.size() >= n_vals ){
				break;
			}else{
				out_vec.push_back(*reinterpret_cast<vbin_t*>((char*)buffer+sizeof(vbin_t)*i));
			}
		}
	}
	
	free(buffer);

	return out_vec;
}


vcov_data::vcov_data(std::string pf, std::string reg){
	open(pf, reg);
}

bool fileExists(const std::string& filepath) {
    if (FILE *file = fopen(filepath.c_str(), "r")) {
        fclose(file);
        return true;
    }
	return false;
}

void vcov_data::open(std::string pf, std::string reg){
	
	file_prefix = pf;
	region = reg;
	
	if( fileExists(file_prefix + ".vcov.bin") ){
		bin.open(file_prefix + ".vcov.bin");
	}else if( fileExists(file_prefix + ".vcov.bin.xz") ){
		bin.open(file_prefix + ".vcov.bin.xz");
	}else if( fileExists(file_prefix + ".vcov.bin.gz") ){
		bin.open(file_prefix + ".vcov.bin.gz");
	}
	
	data_parser dp;
	
	dp.add_field(chr, 0);
	dp.add_field(pos, 1);
	dp.add_field(ref, 2);
	dp.add_field(alt, 3);
	dp.add_field(flipped, 4);
	dp.add_field(mac, 5);
	dp.add_field(var, 6);
	dp.add_field(rid, 7);
	dp.add_field(bin_start_idx, 8);
	dp.add_field(bin_nvals, 9);
	dp.add_matrix(GtU, false, 10);
	
	// dp.parse_file(file_prefix + ".vcov.idx.gz", reg);

	int n_rows = 0;

	indexed_hts_file htsf(file_prefix + ".vcov.idx.gz", reg);
	dp.parse_file(htsf, n_rows);
	htsf.close();
	
}

Eigen::MatrixXd vcov_data::getGtG(int i_s, int i_e){
	int sz = i_e - i_s;
	
	Eigen::MatrixXd out = Eigen::MatrixXd::Zero(sz, sz);
	
	int i = 0; 
	int j = 0;
	
	for(int k = i_s; k < i_e; k++){

		double& mac_k = mac[k];
		long nv = sz < bin_nvals[k] ? sz : bin_nvals[k];
		
		if(  mac_k > 0 ){
			std::vector<vbin_t> vals = bin.get_data(bin_start_idx[k] * sizeof(vbin_t), nv);
			
			j = i;
			for(const vbin_t& v: vals){
				if( j < sz ){
					// out(i,j) = (((double) v) * 2.0 * mac_k)/VBIN_T_MAX;
					out(i,j) = unpack_dp((double) v, mac_k, mac[k + j-i]);
					out(j,i) = out(i,j);
				}else{
					break;
				}
				j++;
			}
		}
		i++;
	}
	
	return out;
}

Eigen::MatrixXd vcov_data::getGtG(const std::vector<int>& v_i, const std::vector<int>& v_j){
	
	int n_row = v_i.size(), n_col = v_j.size();
	int first = v_i[0], last = v_i[0];
	for( const auto& x : v_i){
		if( x < last ){
			std::cerr << "Fatal: Unsorted input in getV.\n";
			abort();
		}
		last = x;
	}
	
	Eigen::MatrixXd out;
	
	if( n_col == 0 ){
		
		out = Eigen::MatrixXd::Zero(n_row, n_row);
		
		int i = 0;
		
		for(const int& ii : v_i ){

			double& mac_ii = mac[ii];
			long nv = v_i[n_row - 1] - v_i[i] + 1 < bin_nvals[ii] ? v_i[n_row - 1] - v_i[i] + 1 : bin_nvals[ii];
			
			if(  mac_ii > 0 ){
				std::vector<vbin_t> vals = bin.get_data(bin_start_idx[ii] * sizeof(vbin_t), nv);
				
				for(int j = i; j < n_row; ++j){
					int jj = v_i[j] - ii;
					if( jj < vals.size() ){
						// out(i,j) = (((double) v) * 2.0 * mac_k)/VBIN_T_MAX;
						out(i,j) = unpack_dp((double) vals[jj], mac_ii, mac[jj+ii]);
						out(j,i) = out(i,j);
					}else{
						break;
					}
				}
			}
			i++;
		}

	}else{
		out = Eigen::MatrixXd::Zero(n_row, n_col);
		
		int i = 0;
		
		for(const int& ii : v_i ){
			out.row(i) = getColGtG(ii, v_j).transpose();
			i++;
		}
	}
	
	return out;
}

Eigen::MatrixXd vcov_data::getV(const std::vector<int>& v_i, const std::vector<int>& v_j){
	
	int n_row = v_i.size(), n_col = v_j.size();
	
	int first = v_i[0], last = v_i[0];
	
	Eigen::MatrixXd out;
	
	bool is_sorted = true;
	
	for( const auto& x : v_i){
		if( x < last ){
			// std::cerr << "Fatal: Unsorted input in getV.\n";
			// abort();
			is_sorted = false;
			// break;
		}
		last = x;
	}
	int nv_total = last - first;
	
	if( n_col == 0 ){
		
		if( n_row >= 250 && n_row >= 0.3*nv_total && is_sorted ){
			
			int s_first = bin_start_idx[first];
			int s_last = bin_start_idx[last];
			
			// std::cerr << "\nLoading vcov data ... \n";
			
			std::vector<vbin_t> vals = bin.get_data(s_first * sizeof(vbin_t), s_last - s_first + 1);
			
			// std::cerr << "Loaded vcov data. Processing matrix ...\n\n";

			out = Eigen::MatrixXd::Zero(n_row, n_row);
	
			int i = 0;
			for(const int& ii : v_i ){

				double& mac_ii = mac[ii];
				long nv = v_i[n_row - 1] - v_i[i] + 1 < bin_nvals[ii] ? v_i[n_row - 1] - v_i[i] + 1 : bin_nvals[ii];
				
				if(  mac_ii > 0 ){
					
					
					for(int j = i; j < n_row; j++){
						int jj = v_i[j] - ii;
						int jv = v_i[j] - ii + (bin_start_idx[ii] - s_first);
						if( ii == ii + jj ){
							out(i,j) = var[ii];
						}else if( jj < vals.size() ){
							// out(i,j) = (((double) v) * 2.0 * mac_k)/VBIN_T_MAX;
							out(i,j) = unpack_dp((double) vals[jv], mac_ii, mac[jj+ii]) - GtU.row(ii).dot(GtU.row(ii+jj));
							out(j,i) = out(i,j);
						}else{
							break;
						}
					}
				}
				i++;
				// if( i % 1000 == 0 ) std::cerr << i << " rows ...\n";
			}
			
			// std::cerr << "Done processing matrix.\n\n";

		}else{
			
			out = Eigen::MatrixXd::Zero(n_row, n_row);
	
			if ( is_sorted ){
				int i = 0;
				for(const int& ii : v_i ){

					double& mac_ii = mac[ii];
					long nv = v_i[n_row - 1] - v_i[i] + 1 < bin_nvals[ii] ? v_i[n_row - 1] - v_i[i] + 1 : bin_nvals[ii];
					
					if(  mac_ii > 0 ){
						std::vector<vbin_t> vals = bin.get_data(bin_start_idx[ii] * sizeof(vbin_t), nv);
						
						for(int j = i; j < n_row; j++){
							int jj = v_i[j] - ii;
							if( ii == ii + jj ){
								out(i,j) = var[ii];
							}else if( jj < vals.size() ){
								// out(i,j) = (((double) v) * 2.0 * mac_k)/VBIN_T_MAX;
								out(i,j) = unpack_dp((double) vals[jj], mac_ii, mac[jj+ii]) - GtU.row(ii).dot(GtU.row(ii+jj));
								out(j,i) = out(i,j);
							}else{
								break;
							}
						}
					}
					i++;
					// if( i % 1000 == 0 ) std::cerr << i << " rows ...\n";
				}

				// std::cerr << "Done processing matrix.\n\n";

			}else{
				int i = 0, j = 0;
				for(int i = 0; i < n_row; i++ ){
					const int& ii = v_i[i];
					out(i,i) = var[ii];
					for( int j = i+1; j < n_row; j++ ){
						const int& jj = v_i[j];
						out(i,j) = getPairV(ii,jj);
						out(j,i) = out(i,j); 
					}
				}
			}
		}
	}else{
		out = Eigen::MatrixXd::Zero(n_row, n_col);
		
		int i = 0;
		
		for(const int& ii : v_i ){
			out.row(i) = getColV(ii, v_j).transpose();
			i++;
		}
	}
	
	return out;
}


Eigen::MatrixXd vcov_data::getColGtG(int ii, const std::vector<int>& v_k ){
	
	int sz = v_k.size();
	
	Eigen::MatrixXd out = Eigen::MatrixXd::Zero(sz, 1);
	
	if(mac[ii] == 0) return out;
	
	int j = 0;
	
	for( const int& k : v_k ){
		if(mac[k] > 0)
		{
			if( k < ii )
			{
				if( ii-k < bin_nvals[k] )
				{
					out(j,0) = unpack_dp((double) bin.get_data( (bin_start_idx[k] +ii-k) * sizeof(vbin_t), 1)[0], mac[k], mac[ii]);
				}
			}else if(k >= ii)
			{
				if( k-ii < bin_nvals[ii] )
				{
					out(j,0) = unpack_dp((double) bin.get_data( (bin_start_idx[ii] + k-ii) * sizeof(vbin_t), 1)[0], mac[k], mac[ii]);
				}
			}
		}
		j++;
	}
	
	return out;
}

Eigen::MatrixXd vcov_data::getColGtG(int ii, int i_s, int i_e){
	
	int sz = i_e - i_s;
	
	Eigen::MatrixXd out = Eigen::MatrixXd::Zero(sz, 1);
	
	if(mac[ii] == 0) return out;
	
	int j = 0;
	
	for(int k = i_s; k < i_e; k++){
		if(mac[k] > 0){
			if( k < ii ){
				if( ii-k < bin_nvals[k] ){
					out(j,0) = unpack_dp((double) bin.get_data( (bin_start_idx[k] +ii-k) * sizeof(vbin_t), 1)[0], mac[k], mac[ii]);
				}
			}else if(k >= ii){
				if( k-ii < bin_nvals[ii] ){
					out(j,0) = unpack_dp((double) bin.get_data( (bin_start_idx[ii] + k-ii) * sizeof(vbin_t), 1)[0], mac[k], mac[ii]);
				}
			}
		}
		j++;
	}
	
	return out;
}

Eigen::MatrixXd vcov_data::getV(int i_s, int i_e){
	Eigen::MatrixXd out = getGtG(i_s, i_e) - GtU.middleRows(i_s,i_e-i_s) * GtU.middleRows(i_s,i_e-i_s).transpose();
	for( int i = i_s, ii = 0; i < i_e; ++i, ++ii){
		for( int j = i, jj = ii; j < i_e; ++j, ++jj ){
			if( i == j ){
				out(ii,jj) = var[i];
			}else{
				if( j - i >= bin_nvals[i] ){
					out(ii,jj) = 0; 
					out(jj,ii) = 0;
				}else if( flipped[i] != flipped[j] ){
					out(ii,jj) = -1.0 * out(ii,jj);
					out(jj,ii) = out(ii,jj);
				}
			}
		}
	}
	return out;
}

Eigen::MatrixXd vcov_data::getColV(int ii, const std::vector<int>& v_k){

	int sz = v_k.size();
	
	Eigen::MatrixXd out = Eigen::MatrixXd::Zero(sz, 1);
	
	if(mac[ii] == 0) return out;
	
	int j = 0;
	for(const int& k : v_k){
		if(mac[k] > 0){
			if( k < ii ){
				if( ii-k < bin_nvals[k] ){
					out(j,0) = unpack_dp((double) bin.get_data( (bin_start_idx[k] +ii-k) * sizeof(vbin_t), 1)[0], mac[k], mac[ii]) - GtU.row(ii).dot(GtU.row(k));
				}
			}else if(k == ii){
				out(j,0) = var[ii];
			}if(k > ii){
				if( k-ii < bin_nvals[ii] ){
					out(j,0) = unpack_dp((double) bin.get_data( (bin_start_idx[ii] + k-ii) * sizeof(vbin_t), 1)[0], mac[k], mac[ii]) - GtU.row(ii).dot(GtU.row(k));
				}
			}
			if( flipped[ii] != flipped[k] ) out(j,0) = (-1.0)*out(j,0);
		}
		j++;
	}
	
	return out;
}

Eigen::MatrixXd vcov_data::getColV(int ii, int i_s, int i_e){
	/*Eigen::MatrixXd out = getColGtG(ii, i_s, i_e) -  GtU.middleRows(i_s,i_e-i_s) * GtU.row(ii).transpose();
	
	int& flip_ii = flipped[ii];
	for(int k = i_s; k < i_e; k++){
		if( flip_ii != flipped[k] ){
			out(k,0) = -1.0*out(k,0);
		}else if( ii == k ){
			out(k,0) = var[ii];
		}else if(){
			
		}
	}*/
	
	
	int sz = i_e - i_s;
	
	Eigen::MatrixXd out = Eigen::MatrixXd::Zero(sz, 1);
	
	if(mac[ii] == 0) return out;
	
	
	for(int k = i_s, j = 0; k < i_e; k++, j++){
		if(mac[k] > 0){
			if( k < ii ){
				if( ii-k < bin_nvals[k] ){
					out(j,0) = unpack_dp((double) bin.get_data( (bin_start_idx[k] +ii-k) * sizeof(vbin_t), 1)[0], mac[k], mac[ii]) - GtU.row(ii).dot(GtU.row(k));
				}
			}else if(k == ii){
				out(j,0) = var[ii];
			}if(k > ii){
				if( k-ii < bin_nvals[ii] ){
					out(j,0) = unpack_dp((double) bin.get_data( (bin_start_idx[ii] + k-ii) * sizeof(vbin_t), 1)[0], mac[k], mac[ii]) - GtU.row(ii).dot(GtU.row(k));
				}
			}
			if( flipped[ii] != flipped[k] ) out(j,0) = (-1.0)*out(j,0);
		}
	}
	
	return out;
}

double vcov_data::getPairV(int ii, int jj){

	if( ii > jj ){
		return getPairV(jj,ii);
	}else if( ii == jj ){
		return var[ii];
	}

	if(mac[ii] == 0) return 0.0;
	if(mac[jj] == 0) return 0.0;
	
	if( jj - ii < bin_nvals[ii] ){
		double out = unpack_dp((double) bin.get_data( (bin_start_idx[ii] + jj-ii) * sizeof(vbin_t), 1)[0], mac[jj], mac[ii]) - GtU.row(ii).dot(GtU.row(jj));
		if( flipped[ii] != flipped[jj] ){
			return (-1)*out;
		}else{
			return out;
		}
	}else{
		return 0.0;
	}
}


Eigen::MatrixXd vcov_data::getV_i(const std::vector<int>& idx, const int offset){
	
	Eigen::MatrixXd out(idx.size(), idx.size());
	
	for( int i = 0; i < idx.size(); i++){
		
		out(i,i) = var[idx[i] + offset];
		
		for( int j = i + 1; j < idx.size(); j++){
			out(i,j) = getPairV(idx[i] + offset, idx[j] + offset);
			out(j,i) = out(i,j);
		}
		
	}
	
	return out;
}

void vcov_data::close(){
	bin.close();
}


void write_vcov_files(genotype_data& g_data, const table& c_data){
	
	const Eigen::SparseMatrix<double> &G = g_data.genotypes;
	
	const Eigen::MatrixXd &X = c_data.data_matrix;
	
	Eigen::MatrixXd U = get_half_hat_matrix(c_data.data_matrix);
	
	Eigen::MatrixXd UtG = (get_half_hat_matrix(c_data.data_matrix).transpose() * G).eval();
	
	//Eigen::MatrixXd& UtG = g_data.UtG;
	
	std::string bin_file_path = global_opts::out_prefix + ".vcov.bin";
	std::string idx_file_path = global_opts::out_prefix + ".vcov.idx.gz";
	
	vcov_bin_gz_out ld_bin_file(bin_file_path);
	BGZF* ld_idx_file = bgzf_open(idx_file_path.c_str(), "w");

	unsigned long int bytes_written = 0;

    vbin_t n_variants = g_data.ld_index_range.size();
	
	//Eigen::MatrixXd ld_matrix = G.transpose() * G;
	
	std::vector<double> macs(G.cols()); 
	for(int i = 0; i < G.cols(); ++i){
		macs[i] = G.col(i).sum(); 
	}
	
	std::string iter_cerr_suffix = " variants ... ";
	
	std::cerr << "Processed LD for ";
	
	int last = 0;
	print_iter_cerr(1, 0, iter_cerr_suffix);
	
	std::string header_string = "# NS=" + std::to_string(G.rows()) + ",NC=" + std::to_string(UtG.rows()) + "\n";
	
	write_to_bgzf(header_string, ld_idx_file);
	
	int n_var = 0;
	for( const auto ld_i : g_data.ld_index_range ){
	
		std::stringstream ld_idx_line;
	
		ld_idx_line.precision(4);
	
		// std::cout << ld_i.first << ", " << ld_i.second << "\n";
	
		double& mac_i = macs[ld_i.first];
		
		int nn = ld_i.second - ld_i.first + 1;
		nn = nn + ld_i.first < G.cols() ? nn : G.cols() - ld_i.first;
		
		Eigen::MatrixXd ld_record;
		
		double var_i = 0;
		
		if( mac_i > 0 ){
			ld_record = ( G.col(ld_i.first).transpose() * G.middleCols(ld_i.first, nn));
			//ld_record = ld_matrix.middleCols(ld_i.first, nn).row(ld_i.first);
			var_i = ld_record(0,0) - UtG.col(ld_i.first).squaredNorm();
		}
		
		g_data.var.push_back(var_i);
		
		ld_idx_line << 
			clean_chrom(g_data.chr[ld_i.first]) << "\t" <<
			g_data.pos[ld_i.first] << "\t" << 
			g_data.ref[ld_i.first] << "\t" << 
			g_data.alt[ld_i.first] << "\t" << 
			g_data.flipped[ld_i.first] << "\t" << 
			mac_i << "\t" << 
			var_i << "\t" << 
			n_var << "\t" <<
			bytes_written/sizeof(vbin_t) << "\t" << nn;
		
		for( int i = 0; i < UtG.rows(); ++i ){
			ld_idx_line << "\t" << UtG(i, ld_i.first);
		}
		ld_idx_line << "\n";
		
		write_to_bgzf(ld_idx_line.str(), ld_idx_file);
		
		for(int i = 0; i < ld_record.cols(); ++i )
		{
			
			// vbin_t tmp_val = (vbin_t) ((unsigned int) round(0.5 * (ld_record(0,i) * VBIN_T_MAX/ mac_i))); 
			
			// vbin_t tmp_val = (vbin_t) ( (unsigned int) ld_record(0,i) );
			vbin_t tmp_val = (vbin_t) pack_dp(ld_record(0,i), mac_i, macs[ld_i.first + i] );

			ld_bin_file.add_to_queue(tmp_val);
			
			bytes_written += sizeof(vbin_t);
		}
		
		n_var++;
		
		// if( n_var % 250 == 0 ) print_iter_cerr(n_var - 250, n_var, iter_cerr_suffix);
		thinned_iter_cerr(last, n_var, iter_cerr_suffix, 250);
	}
	std::cerr << "\rFinished processing LD for " << n_var << " variants.\n";
	
	bgzf_close(ld_idx_file);
	ld_bin_file.close();

	if ( tbx_index_build(idx_file_path.c_str(), 14, &tbx_conf_vcf)!=0 )
	{
		std::cerr << "Fatal error: Couldn't create index " << idx_file_path << ".gz.csi\n";
		exit(1);
	}
	
}
