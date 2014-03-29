#ifndef __HAMMER_IT_CONSENSUS_HPP__
#define __HAMMER_IT_CONSENSUS_HPP__

#include "HSeq.hpp"

#include <limits>
#include <boost/numeric/ublas/matrix.hpp>

namespace hammer {
namespace iontorrent {

inline
std::pair<hammer::HomopolymerRun, double> consensus(const boost::numeric::ublas::matrix<double>& scores) {
  double inf = -std::numeric_limits<double>::infinity();

  double max = inf;
  uint8_t nucl = 0;
  uint8_t len = 1;
  for (uint8_t j = 0; j < 4; ++j)
    for (uint8_t k = 1; k < 64; ++k)
      if (scores(j, k) > max) {
        nucl = j;
        len = k;
        max = scores(j, k);
      }

  return std::make_pair(hammer::HomopolymerRun(nucl, len), max);
}

};
};

#endif // __HAMMER_IT_CONSENSUS_HPP__
