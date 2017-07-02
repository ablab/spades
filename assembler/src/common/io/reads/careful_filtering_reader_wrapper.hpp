//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************
#pragma once
//todo rename file
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
    std::pair<size_t, size_t> c1 = LongestValidCoords(r.first());
    std::pair<size_t, size_t> c2 = LongestValidCoords(r.second());
    size_t len1 = c1.second - c1.first;
    size_t len2 = c2.second - c2.first;
    //FIXME is it the desired logic?
    if (len1 == 0 || len2 == 0) {
        return PairedRead();
    }

    if (len1 == r.first().size() && len2 == r.second().size()) {
        return r;
    }

    return PairedRead(r.first().Substr(c1.first, c1.second),
                      r.second().Substr(c2.first, c2.second),
                      r.orig_insert_size());
}


//todo rewrite without eof
template<typename ReadType>
class CarefulFilteringWrapper : public DelegatingWrapper<ReadType> {
    typedef DelegatingWrapper<ReadType> base;
public:
  /*
   * Default constructor.
   *
   * @param reader Reference to any other reader (child of IReader).
   */
    CarefulFilteringWrapper(typename base::ReadStreamPtrT reader_ptr) :
                base(reader_ptr),
                eof_(false) {
            StepForward();
        }

    bool eof() override {
        return eof_;
    }

  /*
   * Read SingleRead from stream.
   *
   * @param read The SingleRead that will store read * data.
   *
   * @return Reference to this stream.
   */
    CarefulFilteringWrapper& operator>>(ReadType& read) override {
        read = next_read_;
        StepForward();
        return *this;
    }

    void reset() override {
        base::reset();
        eof_ = false;
        StepForward();
    }

private:
    bool eof_;
    ReadType next_read_;

  /*
   * Read next valid read in the stream.
   */
    void StepForward() {
        while (!base::eof()) {
            base::operator >>(next_read_);
            next_read_ = LongestValid(next_read_);
            if (next_read_.IsValid()) {
                return;
            }
        }
        eof_ = true;
    }
};

//FIXME Rename and move!
template<class ReadType>
std::shared_ptr<ReadStream<ReadType>> CarefulFilteringWrap(std::shared_ptr<ReadStream<ReadType>> reader_ptr) {
    //return reader_ptr = make_shared<CarefulFilteringWrapper<ReadType>>(reader_ptr, false, LibraryOrientation::Undefined);
    return std::make_shared<CarefulFilteringWrapper<ReadType>>(reader_ptr);
}

//FIXME Rename and move!
template<class ReadType>
ReadStreamList<ReadType> CarefulFilteringWrap(const ReadStreamList<ReadType>& readers) {
    ReadStreamList<ReadType> answer;
    for (size_t i = 0; i < readers.size(); ++i) {
        answer.push_back(CarefulFilteringWrap<ReadType>(readers.ptr_at(i)));
    }
    return answer;
}

}
