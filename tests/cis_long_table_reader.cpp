#include "cis_long_table_reader.h"
using namespace std;

void CisLongTableReader::load(const string &filepath) {
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
  auto regex_header = regex("#chrom\tpos\tref\talt\tgene\tbeta\tse\tpval");

  bool header_done = false;
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
      auto rec = make_shared<CisLongGeneRecord>();
      rec->chrom = tokens.at(0);
      rec->pos = stoull(tokens.at(1));
      rec->ref = tokens.at(2);
      rec->alt = tokens.at(3);
      rec->gene = tokens.at(4);
      rec->beta = extract_fp<long double>(tokens.at(5));
      rec->se = extract_fp<long double>(tokens.at(6));
      rec->pval = extract_fp<long double>(tokens.at(7));

      // Keys
//      string chrpos = rec->chrom + ":" + to_string(rec->pos);
//      string variant = rec->chrom + ":" + to_string(rec->pos) + "_" + rec->ref + "/" + rec->alt;

      // Insert
      records.push_back(rec);
      this->index.emplace(rec->gene, rec);
    }
  }

  if (records.empty()) {
    throw std::runtime_error("No results were read from file: " + filepath);
  }
}

CisLongTableReader::CisLongTableReader(const string &file) {
  load(file);
}

shared_ptr<CisLongGeneRecord> CisLongTableReader::get_record(const string &i) {
  auto it = this->index.find(i);
  if (it != this->index.end()) {
    return it->second;
  }
  else {
    auto null = shared_ptr<CisLongGeneRecord>();
    return null;
  }
}

bool CisLongTableReader::operator==(const CisLongTableReader& other) {
  if (records.size() != other.records.size()) {
    return false;
  }

  const size_t nrecords = records.size();
  for (size_t i = 0; i < nrecords; i++) {
    auto& rec_a = records[i];
    auto& rec_b = other.records[i];

    if (rec_a->chrom != rec_b->chrom)
      return false;
    if (rec_a->gene != rec_b->gene)
      return false;
    if (rec_a->pos != rec_b->pos)
      return false;
    if (rec_a->ref != rec_b->ref)
      return false;
    if (rec_a->alt != rec_b->alt)
      return false;

    if (!approx_equal(rec_a->beta, rec_b->beta))
      return false;
    if (!approx_equal(rec_a->se, rec_b->se))
      return false;
    if (!approx_equal(rec_a->pval, rec_b->pval, PVALUE_MAX_REL_DIFF))
      return false;
  }

  return true;
}