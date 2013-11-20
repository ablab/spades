#include <iostream>
#include <algorithm>
#include <map>
#include <string>
#include <vector>

#include "utils.hpp"

namespace cclean {

  std::string reverseComplement(const std::string& read) {
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
          score += additional::BruteMatchScore;
        break;
      }
      case 1: { //insert
        for (int i = 0; i < num; ++i, ++query_pos)
          score -= qual[query_pos] / additional::BruteMismatchScore;
        break;
      }
      case 2: { //del
        for (int i = 0; i < num; ++i, ++ref_pos)
          score -= qual[query_pos] / additional::BruteMismatchScore;
        break;
      }
      default:
        break;
    }
  }

  return score;
}

  // end of namespace cclean
}
