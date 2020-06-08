//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************
#pragma once

#include "filtering_reader_wrapper.hpp"
#include "single_read.hpp"
#include "paired_read.hpp"

namespace io {

inline std::pair<size_t, size_t> LongestValidCoords(const SingleRead& r) {
    const size_t none = -1ul;

    size_t best_len = 0;
    size_t best_pos = none;
    size_t pos = none;
    const std::string &seq = r.GetSequenceString();
    size_t sz = seq.length();
    for (size_t i = 0; i <= sz; ++i) {
        if (i < sz && is_nucl(seq[i])) {
            if (pos == none)
                pos = i;
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

inline void LongestValid(SingleRead& r) {
    size_t from, to;
    std::tie(from, to) = LongestValidCoords(r);
    size_t newsize = to - from;
    if (newsize == 0)
        r = SingleRead();
    else if (newsize < r.size())
        r = r.Substr(from, to, /* validate */ false); // FIXME: do inplace
}

inline void LongestValid(PairedRead& r) {
    LongestValid(r.first());
    LongestValid(r.second());
}

template<typename ReadType>
class LongestValidRetainingWrapper : public DelegatingWrapper<ReadType> {
    typedef DelegatingWrapper<ReadType> base;
public:
    LongestValidRetainingWrapper(typename base::ReadStreamT reader_ptr) :
            base(std::move(reader_ptr)) {}

    LongestValidRetainingWrapper& operator>>(ReadType& read) {
        this->reader() >> read;
        LongestValid(read);
        return *this;
    }
};

template<class ReadType>
ReadStream<ReadType> LongestValidWrap(ReadStream<ReadType> reader_ptr) {
    return LongestValidRetainingWrapper<ReadType>(std::move(reader_ptr));
}

template<class ReadType>
ReadStreamList<ReadType> LongestValidWrap(ReadStreamList<ReadType> readers) {
    ReadStreamList<ReadType> answer;
    for (auto &reader : readers)
        answer.push_back(LongestValidWrap<ReadType>(reader));

    return answer;
}

}
