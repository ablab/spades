#ifndef UTILS_H_
#define UTILS_H_

#include <ssw/ssw_cpp.h> // Striped Smith-Waterman aligner
#include <io/read.hpp>
#include "additional.cpp"

namespace cclean {

std::string reverseComplement(const std::string& read);
double GetScoreWithQuality(const StripedSmithWaterman::Alignment &a,
                                            const Quality &qual);
}
#endif /* UTILS_H_ */
