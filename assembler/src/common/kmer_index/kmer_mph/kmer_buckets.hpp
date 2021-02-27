//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "adt/lemiere_mod_reduce.hpp"
#include <iosfwd>
#include <cstdlib>


namespace kmer {

class KMerSegmentPolicyBase {
public:
    explicit KMerSegmentPolicyBase(size_t num_segments = 0)
            : num_segments_(num_segments) {}

    void reset(size_t num_segments) {
        num_segments_ = num_segments;
    }

    size_t num_segments() const { return num_segments_; }

    void BinRead(std::istream& is);
    void BinWrite(std::ostream& os) const;

protected:
    uint64_t reduce(uint64_t val) const {
        return mod_reduce::multiply_high_u64(val, num_segments_);
    }

    size_t num_segments_ = 0;
};


template<class Seq>
class KMerSegmentPolicy : public KMerSegmentPolicyBase {
    typedef typename Seq::hash hash;

  public:
    using KMerSegmentPolicyBase::KMerSegmentPolicyBase;

    size_t operator()(const Seq &s) const {
        if (this->num_segments_ == 1)
            return 0;

        return reduce(hash()(s));
    }

    template<class Ref>
    size_t operator()(Ref s) const {
        if (this->num_segments_ == 1)
            return 0;

        return reduce(hash()(s.data(), s.size()));
    }
};

}
