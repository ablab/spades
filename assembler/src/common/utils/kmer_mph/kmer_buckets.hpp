//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "adt/lemiere_mod_reduce.hpp"
#include <cstdlib>

namespace kmer {

template<class Seq>
class KMerSegmentPolicy {
    typedef typename Seq::hash hash;

public:
    explicit KMerSegmentPolicy(size_t num_segments = 0)
            : num_segments_(num_segments) {}

    void reset(size_t num_segments) {
        num_segments_ = num_segments;
    }

    size_t num_segments() const { return num_segments_; }

    size_t operator()(const Seq &s) const {
        if (num_segments_ == 1)
            return 0;
        
        return mod_reduce::multiply_high_u64(hash()(s), num_segments_);
    }

    template<class Ref>
    size_t operator()(Ref s) const {
        if (num_segments_ == 1)
            return 0;

        return mod_reduce::multiply_high_u64(hash()(s.data(), s.size()), num_segments_);
    }

private:
    size_t num_segments_ = 0;
};

}
