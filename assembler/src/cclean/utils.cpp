#include <iostream>
#include <algorithm>
#include <map>
#include <string>
#include <vector>

#include "utils.hpp"
#include <ssw/ssw_cpp.h>
#include <ssw/ssw_cpp.h> // Striped Smith-Waterman aligner
#include <io/read.hpp>
#include "additional.cpp"

namespace cclean_utils {

inline std::string ReverseComplement(const std::string& read) {
  std::map<char, char> reverse;
  reverse['C'] = 'G';
  reverse['G'] = 'C';
  reverse['T'] = 'A';
  reverse['A'] = 'T';
  reverse['N'] = 'N';

  std::vector<char> res;
  for(int i = 0; i < (int) read.length(); ++i) {
   res.push_back(reverse[read[i]]);
  }

  std::reverse(res.begin(), res.end());
  return std::string(res.begin(), res.end());
}

double GetScoreWithQuality(const StripedSmithWaterman::Alignment &a,
                                            const Quality &qual)
{ // Try to get more realistic align score depend on read quality
  // Mathes and mismatches get from cigar alignment string below
  double score = 0.0;
  int ref_pos = 0, query_pos = 0;
  for (std::vector<uint32_t>::const_iterator it = a.cigar.begin();
       it != a.cigar.end(); ++it) {

    int num = (*it & 0xFFFFFFF0) >> 4;
    int op_code = *it & 0x0000000F;

    switch (op_code) {
      case 0: { //match
        for (int i = 0; i < num; ++i, ++ref_pos, ++query_pos)
          score += additional::MatchScore;
        break;
      }
      case 1: { //insert
        for (int i = 0; i < num; ++i, ++query_pos)
          score -= (double)qual[query_pos] / additional::MismatchScore;
        break;
      }
      case 2: { //del
        for (int i = 0; i < num; ++i, ++ref_pos)
          score -= (double)qual[query_pos] / additional::MismatchScore;
        break;
      }
      default:
        break;
    }
  }
  return score;
}

Read CutRead(const Read &r, int start_pos, int end_pos) {
  if(start_pos > end_pos)  return r;
  //  Step 1: cutting read sequence
  Read read = r;
  std::string read_seq = read.getSequenceString();
  std::string cuted_read_seq(std::string(read_seq, 0, start_pos) +
                             std::string(read_seq, end_pos + 1));
  read.setSequence(cuted_read_seq.c_str());

  //  Step 2: cutting read quality string
  std::string qual_string = read.getQuality().str();
  if(qual_string.empty())  return read;
  std::string cuted_qual_string(std::string(qual_string, 0, start_pos) +
                                std::string(qual_string, end_pos + 1));
  read.setQuality(cuted_qual_string.c_str(), 0);
  return read;
}

void RestoreFromCigar(const std::string& ref, const std::string& query,
                      std::string& out_ref, std::string& out_query,
                      const StripedSmithWaterman::Alignment& a) {

  std::vector<char> aligned_ref, aligned_query;
  int ref_pos = 0, query_pos = 0;
  for (std::vector<uint32_t>::const_iterator it = a.cigar.begin();
       it != a.cigar.end(); ++it) {
    int num = (*it & 0xFFFFFFF0) >> 4;
    int op_code = *it & 0x0000000F;

    switch (op_code) {
      case 0: { //match
        for (int i = 0; i < num; ++i) {
          aligned_ref.push_back(ref[a.ref_begin + ref_pos++]);
          aligned_query.push_back(query[a.query_begin + query_pos++]);
        }
        break;
      }
      case 1: { //insert
        for (int i = 0; i < num; ++i) {
          aligned_ref.push_back('-');
          aligned_query.push_back(query[a.query_begin + query_pos++]);
        }
        break;
      }
      case 2: { //del
        for (int i = 0; i < num; ++i) {
          aligned_ref.push_back(ref[a.ref_begin + ref_pos++]);
          aligned_query.push_back('-');
        }
        break;
     }
      default:
        break;
    }

  }

  out_ref = std::string(aligned_ref.begin(), aligned_ref.end());
  out_query = std::string(aligned_query.begin(), aligned_query.end());
}

std::unordered_map<std::string, std::string> ProcessArgs(int argc, char *argv[],
                                                         bool *ok, std::string *error)
{
  std::unordered_map<std::string, std::string> options;

  // Process arg
  for(int i = 1; i < argc; ++i)
  {
    std::string arg;
    int j = 0;
    char next_ch = argv[i][j];

    while (next_ch != ':' && next_ch != '=') {
      if (next_ch == '\0') {
        (*ok) = false;
        (*error) = "Bad argument " + std::string(argv[i]);
        return options;
      }
      arg += next_ch;
      ++j;
      next_ch = argv[i][j];
    }

    // Process argument value
    ++j;
    next_ch = argv[i][j];

    std::string val;
    while (next_ch != '\0') {
      val += next_ch;
      ++j;
      next_ch = argv[i][j];
    }

    if (arg.empty() || val.empty()) {
      if (next_ch == '\0') {
        (*ok) = false;
        (*error) = "Bad argument " + std::string(argv[i]);
        return options;
      }
    }
    options[arg] = val;
  }

  // Process sinonims
  std::unordered_map<std::string, std::vector<std::string> > reqParams;
  // Add new required args here with sinonims
  reqParams["mode"] = {"--m", "--mode", "--MODE"};
  reqParams["config"] = {"--c", "--config", "--CONFIG"};
  reqParams["input"] = {"--i", "--input", "--INPUT"};
  reqParams["output"] = {"--o", "--output", "--OUTPUT"};
  reqParams["database"] = {"--d", "--database", "--DATABASE"};
  // not required args here
  std::unordered_map<std::string, std::vector<std::string> > notReqParams;
  notReqParams["mlen"] = {"--ml", "--mlen", "--MLEN"};
  notReqParams["inform"] = {"--in", "--inform", "--INFORM"};

  for (auto kv: options) {
    std::string arg = kv.first;
    bool correct_arg = false;
    for (auto sinonim: reqParams) {
      if (std::find(sinonim.second.begin(), sinonim.second.end(), arg) !=
          sinonim.second.end()) {
        options.erase(kv.first);
        options[sinonim.first] = kv.second;
        correct_arg = true;
        break;
      }
      if (sinonim.first == kv.first)
        correct_arg = true;
    }

    for (auto sinonim: notReqParams) {
      if (std::find(sinonim.second.begin(), sinonim.second.end(), arg) !=
          sinonim.second.end()) {
        options.erase(kv.first);
        options[sinonim.first] = kv.second;
        correct_arg = true;
        break;
      }
      if (sinonim.first == kv.first)
        correct_arg = true;
    }
    if (!correct_arg) {
      (*ok) = false;
      (*error) = "Bad argument " + std::string(kv.first);
      return options;
    }
  }
  // Check that all required args was specified
  for (auto reqArg: reqParams) {
    std::string arg = reqArg.first;
    bool hasArg = false;
    for (auto opt: options)
      if (opt.first == arg) {
        hasArg = true;
        break;
      }
    if (!hasArg) {
      (*ok) = false;
      (*error) = "Missing argument " + arg;
      return options;
    }
  }

  (*ok) = true;
  return options;
}

  // end of namespace cclean_utils
}
