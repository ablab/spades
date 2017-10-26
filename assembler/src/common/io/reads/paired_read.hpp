//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "single_read.hpp"

#include <string>
#include <utility>

namespace io {

template<class SingleReadType>
class UniversalPairedRead {
public:
    typedef SingleReadType SingleReadT;
    typedef int32_t InsertSizeT;

    //Prefer using creating new read with constructor to the next three methods
    SingleReadT &first() {
        return first_;
    }

    SingleReadT &second() {
        return second_;
    }

    void set_orig_insert_size(size_t is) {
        insert_size_ = is;
    }

    const SingleReadT &first() const {
        return first_;
    }

    const SingleReadT &second() const {
        return second_;
    }

    size_t orig_insert_size() const {
        return insert_size_;
    }

    InsertSizeT corrected_insert_size() const {
        return InsertSizeT(insert_size_ - first_.GetLeftOffset() - second_.GetRightOffset());
    }

    InsertSizeT distance() const {
        return corrected_insert_size() - (InsertSizeT) second_.size();
    }

    InsertSizeT gap() const {
        return corrected_insert_size() - (InsertSizeT) nucl_count();
    }

    size_t size() const {
        return std::max(first_.size(), second_.size());
    }

    size_t nucl_count() const {
        return first_.size() + second_.size();
    }

    bool IsValid() const {
        return first_.IsValid() && second_.IsValid();
    }

    const SingleReadT &operator[](size_t i) const {
        if (i == 0) {
            return first_;
        } else if (i == 1) {
            return second_;
        }
        VERIFY(false);
        return first_;
    }

    const UniversalPairedRead operator!() const {
        return UniversalPairedRead(!second_, !first_, insert_size_);
    }

    bool operator==(const UniversalPairedRead &paired_read) const {
        return first_ == paired_read.first_ &&
               second_ == paired_read.second_ &&
               insert_size_ == paired_read.insert_size_;
    }

    bool BinRead(std::istream &file, size_t estimated_is) {
        first_.BinRead(file);
        second_.BinRead(file);

        insert_size_ = estimated_is;
        return !file.fail();
    }

    bool BinWrite(std::ostream &file, bool rc1 = false, bool rc2 = false) const {
        first_.BinWrite(file, rc1);
        second_.BinWrite(file, rc2);

        return !file.fail();
    }

    UniversalPairedRead() : first_(), second_(), insert_size_(0) { }

    UniversalPairedRead(const SingleReadT &first,
                        const SingleReadT &second,
                        size_t insert_size)
            : first_(first), second_(second), insert_size_(insert_size) { }


private:
    SingleReadT first_;
    SingleReadT second_;
    size_t insert_size_;
};

typedef UniversalPairedRead<SingleRead> PairedRead;

inline std::ostream &operator<<(std::ostream &os, const PairedRead &read) {
    os << "Single read first=" << read.first() << " second=" << read.second() << std::endl;
    return os;
}

typedef UniversalPairedRead<SingleReadSeq> PairedReadSeq;

inline std::ostream &operator<<(std::ostream &os, const PairedReadSeq &read) {
    os << "Paired read first=" << read.first() << " second=" << read.second() << std::endl;
    return os;
}

}
