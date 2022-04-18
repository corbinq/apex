#ifndef APEX_META_SVAR_READER_H
#define APEX_META_SVAR_READER_H

#include <string>
#include <map>
#include <vector>
#include <memory>
#include "reader_util.h"

struct SvarRecord {
  std::string chrom;
  uint64_t pos;
  std::string ref;
  std::string alt;
  std::string gene;
  std::vector<std::string> studies;
  double beta;
  double se;
  long double pval;
  double log_pval;
};

/**
 * Class to read single variant meta-analysis files.
 */
class SvarReader {
protected:
  std::vector<std::shared_ptr<SvarRecord>> records;
  std::map<std::string, std::shared_ptr<SvarRecord>> index;
public:
  SvarReader(const std::string &file);
  void load(const std::string &file);
  std::shared_ptr<SvarRecord> get_record(const std::string &i);

  bool operator==(const SvarReader& other);
};

#endif //APEX_META_SVAR_READER_H