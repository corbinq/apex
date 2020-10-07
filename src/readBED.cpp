#include "readBED.hpp"


void bed_data::readBedHeader(const char *file_name)
{
	if( ids.file.size() != 0 ) return;
	
	htsFile *htsf = hts_open(file_name, "r");
	
	kstring_t str = {0,0,0};

	int line_no = 0;

	while( hts_getline(htsf, KS_SEP_LINE, &str) )
	{
		line_no++;
		if ( str.s[0] == '#' && str.s[1] == '#' ) continue;
		
		header_line = line_no;
		int n_fields;
		int *offsets = ksplit(&str, 0, &n_fields);
		for( int i = 4; i < n_fields; i++)
		{
			ids.file.push_back(std::string(str.s + offsets[i]));
		}
		break;
	}
	
	ks_free(&str);
	hts_close(htsf);

}

void bed_data::readBedFile(const char *file_name, const std::vector<std::string>& regions)
{
	if( ids.file.size() == 0 ){
		readBedHeader(file_name);
	}
	if( ids.keep.size() == 0 ){
		ids.setKeepIDs(ids.file);
	}
	
	std::string file_str = std::string(file_name);
	
	indexed_hts_file htsf(file_str, regions);

	int line_no = 0;
	n_genes = 0;
	
	data_matrix = Eigen::MatrixXd(ids.idx.size(), 50000);

	data_parser dp;
	dp.add_field(chr,0);
	dp.add_field(start,1);
	dp.add_field(end,2);
	dp.add_field(gene_id,3);
	dp.add_matrix(data_matrix, true, ids.idx, 4);
	
	kstring_t str = {0,0,0};

	while( htsf.next_line(str) >= 0 )
	{
		if( !str.l ) break;
		
		line_no++;
		
		if ( str.s[0] == '#' ){
			continue;
		}else if( regions.size() == 0 && line_no <= header_line){
			continue;
		}
		
		dp.parse_fields(str);
		
		if( global_opts::trim_gene_ids )
		{
			gene_id[n_genes] = gene_id[n_genes].substr(0, gene_id[n_genes].find("."));
		}
		
		n_genes++;
		
		if(n_genes % 500  == 0 ) {
			std::cerr << "Processed expression data for " << n_genes << " genes ... \r";
		}
	}
	
	ks_free(&str);
	
	htsf.close();

	data_matrix.conservativeResize(Eigen::NoChange, n_genes);
}

