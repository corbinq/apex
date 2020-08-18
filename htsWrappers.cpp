#include "htsWrappers.hpp"

using namespace std;

void indexed_hts_file::open(const string& prefix, const vector<string>& reg)
{
	regions = reg;
	
	htsf = hts_open(prefix.c_str(), "r");
    tbx = tbx_index_load(prefix.c_str());
	
	itr = NULL;
	
	c_reg = 0;
	n_reg = regions.size();
	
    if( n_reg == 0 ){
		const char** chroms = tbx_seqnames(tbx, &n_reg);
		for( int i = 0; i < n_reg; ++i ){
			regions.push_back(chroms[i]);
		}
		free(chroms);
	}
}

int indexed_hts_file::next_line(kstring_t& str)
{
	if( !itr ){
		itr = tbx_itr_querys(tbx, regions[c_reg].c_str() );
	}
	if( tbx_itr_next(htsf, tbx, itr, &str) < 0 ){
		if( c_reg + 1 < n_reg ){
			c_reg++;
			itr = tbx_itr_querys(tbx, regions[c_reg].c_str() );
			return tbx_itr_next(htsf, tbx, itr, &str);
		}else{
			return -1;
		}
	}
	return 1;
}

void bcf_seek(bcf_srs_t* sr, const string& chr, const int& start, const int& end){
	
	string region = chr;
	
	if( start > 0 ){
		region += ":" + to_string(start);
	}
	if( end > 0 ){
		region += "-" + to_string(end);
	}
	
	bcf_sr_regions_t *reg = bcf_sr_regions_init(region.c_str(),0,0,1,-2);
	
	bcf_sr_seek(sr, reg->seq_names[0], reg->start);
	
	return;
}

void bgzip_file(string file_name, int index){
	
	// index == 0 => Do not build any index. 
	// index == 1 => Build gzi byte index. 
	// index == 2 => Build csi position-based index. 
	
	BGZF *fp;
    void *buffer;
	
	int f_src = open(file_name.c_str(), O_RDONLY);
	
	string file_name_gz = file_name + ".gz";
	fp = bgzf_open(file_name_gz.c_str(), "w\0");

	if( index == 1 ) bgzf_index_build_init(fp);
	
	buffer = malloc(BUFFER_SIZE);

	int c;
	
	while ((c = read(f_src, buffer, BUFFER_SIZE)) > 0){
		if (bgzf_write(fp, buffer, c) < 0){
			cerr << "Fatal error: Couldn't write to " << file_name << ".gz\n";
			exit(1);
		} 
	}

	if( index == 1 ){
		if (bgzf_index_dump(fp, file_name.c_str(), ".gz.gzi") < 0){
			cerr << "Fatal error: Couldn't create index " << file_name << ".gz.gzi\n";
			exit(1);
		}
	}
	
	if (bgzf_close(fp) < 0){
			cerr << "Fatal error: Couldn't close " << file_name << ".gz\n";
			exit(1);
	}
	
	if( index == 2 ){
		if ( tbx_index_build(file_name_gz.c_str(), 14, &tbx_conf_vcf)!=0 )
		{
			cerr << "Fatal error: Couldn't create index " << file_name << ".gz.csi\n";
			exit(1);
		}
	}
	
	unlink(file_name.c_str());
	free(buffer);
	close(f_src);
}

vector<string> get_chroms(string file_name, vector<int>& n_variants){
	
	// get a list of chromosomes in the file
	// this code is directly adapted from bcftools:
	// https://github.com/samtools/bcftools/blob/develop/vcfindex.c
	
	n_variants.clear();
	vector<string> chroms;
	
	const char **seq;
    int nseq, stats;
    tbx_t *tbx = NULL;
    hts_idx_t *idx = NULL;

    htsFile *fp = hts_open(file_name.c_str(),"r");
    if ( !fp ) { fprintf(stderr,"Could not read %s\n", file_name.c_str()); abort(); }
    bcf_hdr_t *hdr = NULL;

	if ( hts_get_format(fp)->format==bcf )
    {
		hdr = bcf_hdr_read(fp);
		if ( !hdr ) { fprintf(stderr,"Could not read the header: %s\n", file_name.c_str()); abort(); }
        idx = bcf_index_load(file_name.c_str());
        if ( !idx ) { fprintf(stderr,"Could not load index for %s\n", file_name.c_str()); abort(); }
    }else{
		if( hts_get_format(fp)->format==vcf ){
			hdr = bcf_hdr_read(fp);
			if ( !hdr ) { fprintf(stderr,"Could not read the header: %s\n", file_name.c_str()); abort(); }
		}
		tbx = tbx_index_load(file_name.c_str());
        if ( !tbx ) { fprintf(stderr,"Could not load index for %s\n", file_name.c_str()); abort(); }
	}
	
	seq = tbx ? tbx_seqnames(tbx, &nseq) : bcf_index_seqnames(idx, hdr, &nseq);
    for (int i=0; i<nseq; i++)
    {
        uint64_t records, v;
        hts_idx_get_stat(tbx ? tbx->idx : idx, i, &records, &v);
        if (stats&2 || !records) continue;
        chroms.push_back(string(seq[i]));
		n_variants.push_back( (int) records );
    }
    free(seq);
	
    hts_close(fp);
    bcf_hdr_destroy(hdr);
    if (tbx) tbx_destroy(tbx);
    if (idx) hts_idx_destroy(idx);
	
    return chroms;
}
