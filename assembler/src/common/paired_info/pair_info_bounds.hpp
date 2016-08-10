//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef OMNI_UTILS_HPP_
#define OMNI_UTILS_HPP_

#include "utils/standard_base.hpp"

namespace omnigraph {


inline size_t PairInfoPathLengthUpperBound(size_t k, size_t insert_size,
                                           double delta) {
    double answer = 0. + (double) insert_size + delta - (double) k - 2.;
    VERIFY(math::gr(answer, 0.));
    return (size_t)std::floor(answer);
}

inline size_t PairInfoPathLengthLowerBound(size_t k, size_t l1, size_t l2,
                                           int gap, double delta) {
    double answer = 0. + (double) gap + (double) k + 2. - (double) l1 - (double) l2 - delta;
    return math::gr(answer, 0.) ? (size_t)std::floor(answer) : 0;
}

}
#endif /* OMNI_UTILS_HPP_ */
