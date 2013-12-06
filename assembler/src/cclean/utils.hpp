#ifndef UTILS_H_
#define UTILS_H_

#include <ssw/ssw_cpp.h> // Striped Smith-Waterman aligner
#include <io/read.hpp>
#include "additional.cpp"

namespace cclean_utils {

std::string ReverseComplement(const std::string& read);

double GetScoreWithQuality(const StripedSmithWaterman::Alignment &a,
                                            const Quality &qual);

// Cut read from start to end position of best aligment with adapter
Read CutRead(const Read &r, int start_pos, int end_pos);

// end of namespace
}
#endif /* UTILS_H_ */
