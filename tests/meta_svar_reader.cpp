#include "meta_svar_reader.h"
#include <fstream>
#include <string>
#include <iostream>
#include <regex>
#include <optional>
#include "reader_util.h"

using namespace std;

void SvarReader::load(const string &file) {
  if (!filepath_exists(file)) {
    throw std::runtime_error("Could not find file: " + file);
  }

  ifstream input_file(file);
  string line;
  auto line_separator = regex("[ \t]");
  auto field_separator = regex(",");

  // Line regexes
  auto regex_header = regex("#chr\tpos\tref\talt\tgene\tstudies\tbeta\tse\tpval");

  bool header_done = false;
  while (getline(input_file, line)) {
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
      auto rec = make_shared<SvarRecord>();
      rec->chrom = tokens.at(0);
      rec->pos = stoull(tokens.at(1));
      rec->ref = tokens.at(2);
      rec->alt = tokens.at(3);
      rec->gene = tokens.at(4);
      copy(sregex_token_iterator(tokens.at(5).begin(), tokens.at(5).end(), field_separator, -1), sregex_token_iterator(), back_inserter(rec->studies));
      rec->beta = extract_fp<double>(tokens.at(6));
      rec->se = extract_fp<double>(tokens.at(7));
      rec->pval = extract_fp<long double>(tokens.at(8));

      if (tokens.size() > 9) {
        rec->log_pval = extract_fp<double>(tokens.at(9));
      }
      else {
        rec->log_pval = 0.0;
      }

      // Keys
//      string chrpos = rec->chrom + ":" + to_string(rec->pos);
//      string variant = rec->chrom + ":" + to_string(rec->pos) + "_" + rec->ref + "/" + rec->alt;

      // Insert
      records.push_back(rec);
      this->index.emplace(rec->gene, rec);
    }
  }

  if (records.empty()) {
    throw std::runtime_error("No single variant meta-analysis results were read from file: " + file);
  }
}

SvarReader::SvarReader(const string &file) {
  load(file);
}

shared_ptr<SvarRecord> SvarReader::get_record(const string &i) {
  auto it = this->index.find(i);
  if (it != this->index.end()) {
    return it->second;
  }
  else {
    auto null = shared_ptr<SvarRecord>();
    return null;
  }
}

bool SvarReader::operator==(const SvarReader& other) {
  if (records.size() != other.records.size()) {
    return false;
  }

  const size_t nrecords = records.size();
  for (size_t i = 0; i < nrecords; i++) {
    auto& rec_a = records[i];
    auto& rec_b = other.records[i];

    if (rec_a->gene != rec_b->gene)
      return false;
    if (rec_a->chrom != rec_b->chrom)
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
    if (!approx_equal(rec_a->log_pval, rec_b->log_pval))
      return false;
  }

  return true;
}