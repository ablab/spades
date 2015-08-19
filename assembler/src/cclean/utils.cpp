//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

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
          score += MatchScore;
        break;
      }
      case 1: { //insert
        for (int i = 0; i < num; ++i, ++query_pos)
          score -= (double)qual[query_pos] / MismatchScore;
        break;
      }
      case 2: { //del
        for (int i = 0; i < num; ++i, ++ref_pos)
          score -= (double)qual[query_pos] / MismatchScore;
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

  // end of namespace cclean_utils
}
