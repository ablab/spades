//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************
#pragma once
//todo rename file
#include "io/reads/delegating_reader_wrapper.hpp"
#include "pipeline/library.hpp"

namespace io {

const size_t none = -1ul;

inline std::pair<size_t, size_t> LongestValidCoords(const SingleRead& r) {
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

inline SingleRead LongestValid(const SingleRead& r,
                               bool /*use_orientation*/ = false,
                               LibraryOrientation /*orientation*/ = LibraryOrientation::FR) {

    std::pair<size_t, size_t> p = LongestValidCoords(r);
    return r.Substr(p.first, p.second);
}

inline PairedRead LongestValid(const PairedRead& r,
                               bool use_orientation = false,
                               LibraryOrientation orientation = LibraryOrientation::FR) {
    std::pair<size_t, size_t> c1 = LongestValidCoords(r.first());
    std::pair<size_t, size_t> c2 = LongestValidCoords(r.second());
    size_t len1 = c1.second - c1.first;
    size_t len2 = c2.second - c2.first;
    if (len1 == 0 || len2 == 0) {
        return PairedRead();
    }
    if (len1 == r.first().size() && len2 == r.second().size()) {
        return r;
    }

    size_t is;
    if (!use_orientation) {
        is = r.insert_size() - c1.first - r.second().size() + c2.second;
    }
    else {
        switch (orientation) {
        case LibraryOrientation::FF:  {
            is = r.insert_size() - c1.first - r.second().size() + c2.second;
            break;
        }
        case LibraryOrientation::RR:  {
            is = r.insert_size() - r.first().size() + c1.second  - c2.first;
            break;
        }
        case LibraryOrientation::FR:  {
            is = r.insert_size() - c1.first - c2.first;
            break;
        }
        case LibraryOrientation::RF:  {
            is = r.insert_size() - r.first().size() + c1.second - r.second().size() + c2.second;
            break;
        }
        default: {
            is = r.insert_size() - c1.first - r.second().size() + c2.second;
            break;
        }
        }
    }

    return PairedRead(r.first().Substr(c1.first, c1.second), r.second().Substr(c2.first, c2.second), is);
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
    CarefulFilteringWrapper(typename base::ReadStreamPtrT reader_ptr,
                                     bool use_orientation = false,
                                     LibraryOrientation orientation = LibraryOrientation::Undefined) :
                base(reader_ptr),
                eof_(false),
                use_orientation_(use_orientation),
                orientation_(orientation) {
            StepForward();
        }

    /* virtual */ bool eof() {
        return eof_;
    }

  /*
   * Read SingleRead from stream.
   *
   * @param read The SingleRead that will store read * data.
   *
   * @return Reference to this stream.
   */
    /* virtual */ CarefulFilteringWrapper& operator>>(ReadType& read) {
        read = next_read_;
        StepForward();
        return *this;
    }

    /* virtual */
    void reset() {
        base::reset();
        eof_ = false;
        StepForward();
    }

private:
    bool eof_;
    bool use_orientation_;
    LibraryOrientation orientation_;
    ReadType next_read_;

  /*
   * Read next valid read in the stream.
   */
    void StepForward() {
        while (!base::eof()) {
            base::operator >>(next_read_);
            next_read_ = LongestValid(next_read_, use_orientation_, orientation_);
            if (next_read_.IsValid()) {
                return;
            }
        }
        eof_ = true;
    }
};

template<class ReadType>
std::shared_ptr<ReadStream<ReadType>> CarefulFilteringWrap(std::shared_ptr<ReadStream<ReadType>> reader_ptr,
                                                           bool use_orientation = false,
                                                           LibraryOrientation orientation = LibraryOrientation::Undefined) {
    //return reader_ptr = make_shared<CarefulFilteringWrapper<ReadType>>(reader_ptr, false, LibraryOrientation::Undefined);
    return std::shared_ptr<CarefulFilteringWrapper<ReadType> >(
            new CarefulFilteringWrapper<ReadType>(reader_ptr, use_orientation, orientation));
}

template<class ReadType>
ReadStreamList<ReadType> CarefulFilteringWrap(const ReadStreamList<ReadType>& readers,
                                              bool use_orientation = false,
                                              LibraryOrientation orientation = LibraryOrientation::Undefined) {
    ReadStreamList<ReadType> answer;
    for (size_t i = 0; i < readers.size(); ++i) {
        answer.push_back(CarefulFilteringWrap<ReadType>(readers.ptr_at(i), use_orientation, orientation));
    }
    return answer;
}

}
