//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************
#pragma once
#include "filtering_reader_wrapper.hpp"
#include "pipeline/library.hpp"

namespace io {

inline std::pair<size_t, size_t> LongestValidCoords(const SingleRead& r) {
    const size_t none = -1ul;

    size_t best_len = 0;
    size_t best_pos = none;
    size_t pos = none;
    std::string seq = r.GetSequenceString();
    for (size_t i = 0; i <= seq.size(); ++i) {
        if (i < seq.size() && is_nucl(seq[i])) {
            if (pos == none) {
                pos = i;
            }
        } else {
            if (pos != none) {
                size_t len = i - pos;
                if (len > best_len) {
                    best_len = len;
                    best_pos = pos;
                }
            }
            pos = none;
        }
    }
    if (best_len == 0) {
        return std::make_pair(0, 0);
    }
    return std::make_pair(best_pos, best_pos + best_len);
}

inline SingleRead LongestValid(const SingleRead& r) {
    std::pair<size_t, size_t> p = LongestValidCoords(r);
    if (p.first == p.second)
        return SingleRead();
    return r.Substr(p.first, p.second);
}

inline PairedRead LongestValid(const PairedRead& r) {
    auto c1 = LongestValidCoords(r.first());
    auto c2 = LongestValidCoords(r.second());

    return PairedRead(r.first().Substr(c1.first, c1.second),
                      r.second().Substr(c2.first, c2.second),
                      r.orig_insert_size());
}


template<typename ReadType>
class LongestValidRetainingWrapper : public DelegatingWrapper<ReadType> {
    typedef DelegatingWrapper<ReadType> base;
public:
    LongestValidRetainingWrapper(typename base::ReadStreamPtrT reader_ptr) :
                base(reader_ptr) {
    }

    LongestValidRetainingWrapper& operator>>(ReadType& read) override {
        this->reader() >> read;
        read = LongestValid(read);
        return *this;
    }
};

template<class ReadType>
std::shared_ptr<ReadStream<ReadType>> LongestValidWrap(std::shared_ptr<ReadStream<ReadType>> reader_ptr) {
    return FilteringWrap<ReadType>(std::make_shared<LongestValidRetainingWrapper<ReadType>>(reader_ptr));
}

template<class ReadType>
ReadStreamList<ReadType> LongestValidWrap(const ReadStreamList<ReadType> &readers) {
    ReadStreamList<ReadType> answer;
    for (size_t i = 0; i < readers.size(); ++i) {
        answer.push_back(LongestValidWrap<ReadType>(readers.ptr_at(i)));
    }
    return FilteringWrap<ReadType>(answer);
}

}
