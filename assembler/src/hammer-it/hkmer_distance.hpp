#ifndef __HAMMER_HKMER_DISTANCE_HPP__
#define __HAMMER_HKMER_DISTANCE_HPP__

#include "hkmer.hpp"
#include "err_helper_table.hpp"

namespace hammer {

static inline size_t hamdistKMer(hammer::HKMer x, hammer::HKMer y,
                                 unsigned tau = -1) {
  unsigned dist = 0;

  for (size_t i = 0; i < hammer::K; ++i) {
    hammer::HomopolymerRun cx = x[i], cy = y[i];
    if (cx.raw != cy.raw) {
      dist += (cx.nucl == cy.nucl ?
               abs(cx.len - cy.len) :  cx.len + cy.len);
      if (dist > tau)
        return dist;
    }
  }

  return dist;
}

static inline size_t distanceHKMer(const hammer::HKMer& x, 
                                   const hammer::HKMer& y,
                                   unsigned tau = -1) {
  unsigned dist = 0;

  size_t i = 0, j = 0, ix = 0, iy = 0;
  hammer::HomopolymerRun cx = x[0], cy = y[0];

  while (i < hammer::K && j < hammer::K) {
    if (ix != i)
      cx = x[ix = i];

    if (iy != j)
      cy = y[iy = j];

    if (cx.raw == cy.raw) {
      ++i, ++j;
      continue;
    }

    if (cx.nucl == cy.nucl) {

      if (cx.len >= 4 && cy.len >= 4) {
        ++i, ++j;
        dist += abs(cx.len - cy.len);
        if (dist > tau)
          return dist;
      } else {
        --cx.len;
        --cy.len;
      }

    } else {

      dist += 1;
      if (dist > tau)
        return dist;

      using namespace hammer::errHelper;
      auto hint = getHint(x, y, i, j, cx.len, cy.len);
      
      switch (hint) {
        case kMismatch:
          --cx.len, --cy.len;
          break;
        case kInsertion:
          --cy.len;
          break;
        case kDeletion:
          --cx.len;
          break;
      }
    }

    if (cx.len == 0) 
      ++i;
    if (cy.len == 0)
      ++j;
  }

  return dist;
}

};
#endif // __HAMMER_HKMER_DISTANCE_HPP__
