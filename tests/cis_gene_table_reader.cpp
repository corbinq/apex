#include "cis_gene_table_reader.h"
using namespace std;

void CisGeneTableReader::load(const string &filepath) {
  if (!filepath_exists(filepath)) {
    throw std::runtime_error("Could not find file: " + filepath);
  }

  unique_ptr<istream> file;
  boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
  ifstream fs(filepath, ios_base::in | ios_base::binary);

  inbuf.push(boost::iostreams::gzip_decompressor());
  inbuf.push(fs);
  file = make_unique<istream>(&inbuf);

  string line;
  auto line_separator = regex("[ \t]");
  auto field_separator = regex(",");

  // Line regexes
  auto regex_header = regex("#chrom\tstart\tend\tgene\tegene_pval.*");

  bool header_done = false;
  bool parse_trait = false;
  while (getline(*file, line)) {
    smatch match;
    if (!header_done) {
      if (regex_search(line, match, regex_header)) {
        header_done = true;
      }
    }
    else {
      // Begin parsing record row
      vector<string> tokens;
      copy(sregex_token_iterator(line.begin(), line.end(), line_separator, -1), sregex_token_iterator(), back_inserter(tokens));

      if (line.substr(0,1) == "#") { continue; }

      // Create record
      auto rec = make_shared<CisGeneRecord>();
      rec->chrom = tokens.at(0);
      rec->start = stoull(tokens.at(1));
      rec->end = stoull(tokens.at(2));
      rec->gene = tokens.at(3);
      rec->egene_pval = extract_fp<long double>(tokens.at(4));
      rec->n_samples = stoull(tokens.at(5));
      rec->n_covar = stoull(tokens.at(6));
      rec->resid_sd = extract_fp<double>(tokens.at(7));
      rec->n_cis_variants = stoull(tokens.at(8));

      // Keys
//      string chrpos = rec->chrom + ":" + to_string(rec->pos);
//      string variant = rec->chrom + ":" + to_string(rec->pos) + "_" + rec->ref + "/" + rec->alt;

      // Insert
      records.push_back(rec);
      this->index.emplace(rec->gene, rec);
    }
  }

  if (records.empty()) {
    throw std::runtime_error("No test results were read from file: " + filepath);
  }
}

CisGeneTableReader::CisGeneTableReader(const string &file) {
  load(file);
}

shared_ptr<CisGeneRecord> CisGeneTableReader::get_record(const string &i) {
  auto it = this->index.find(i);
  if (it != this->index.end()) {
    return it->second;
  }
  else {
    auto null = shared_ptr<CisGeneRecord>();
    return null;
  }
}

bool CisGeneTableReader::operator==(const CisGeneTableReader& other) {
  if (records.size() != other.records.size()) {
    return false;
  }

  const size_t nrecords = records.size();
  for (size_t i = 0; i < nrecords; i++) {
    auto& rec_a = records[i];
    auto& rec_b = other.records[i];

    if (rec_a->gene != rec_b->gene)
      return false;
    if (rec_a->start != rec_b->start)
      return false;
    if (rec_a->end != rec_b->end)
      return false;
    if (rec_a->n_samples != rec_b->n_samples)
      return false;
    if (rec_a->n_covar != rec_b->n_covar)
      return false;
    if (rec_a->n_cis_variants != rec_b->n_cis_variants)
      return false;
    if (!approx_equal(rec_a->egene_pval, rec_b->egene_pval, PVALUE_MAX_REL_DIFF))
      return false;
    if (!approx_equal(rec_a->resid_sd, rec_b->resid_sd))
      return false;
  }

  return true;
}