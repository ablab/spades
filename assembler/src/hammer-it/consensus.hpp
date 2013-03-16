#ifndef __HAMMER_IT_CONSENSUS_HPP__
#define __HAMMER_IT_CONSENSUS_HPP__

#include <limits>
#include <boost/numeric/ublas/matrix.hpp>

namespace hammer {
namespace iontorrent {

std::pair<hammer::HomopolymerRun, double> consensus(const boost::numeric::ublas::matrix<double>& scores) {
  double inf = -std::numeric_limits<double>::infinity();

  unsigned nucl = 0, len = 1; double max = inf;
  for (unsigned j = 0; j < 4; ++j)
    for (unsigned k = 1; k < 64; ++k)
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
