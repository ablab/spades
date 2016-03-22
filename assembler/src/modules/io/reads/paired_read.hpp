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


class PairedRead {
public:
    typedef SingleRead SingleReadT;
    typedef int16_t size_type;

    PairedRead() : first_(), second_(), insert_size_(0) { }

    PairedRead(const SingleRead &first,
               const SingleRead &second,
               size_t insert_size)
            : first_(first), second_(second), insert_size_(insert_size) { }

    const SingleRead &first() const {
        return first_;
    }

    const SingleRead &second() const {
        return second_;
    }

    size_t insert_size() const {
        return insert_size_;
    }

    size_t distance() const {
        return insert_size_ - second_.size();
    }

    size_t gap() const {
        return insert_size_ - first_.size() - second_.size();
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

    const SingleRead &operator[](size_t i) const {
        if (i == 0) {
            return first_;
        } else if (i == 1) {
            return second_;
        }
        VERIFY(false);
        return first_;
    }

    const PairedRead operator!() const {
        return PairedRead(!second_, !first_, insert_size_);
    }

    bool operator==(const PairedRead &pairedread) const {
        return first_ == pairedread.first_ &&
               second_ == pairedread.second_ &&
               insert_size_ == pairedread.insert_size_;
    }

    bool BinWrite(std::ostream &file, bool rc1 = false, bool rc2 = false) const {
        first_.BinWrite(file, rc1);
        second_.BinWrite(file, rc2);

        return !file.fail();
    }

    void print_size() const {
        first_.print_size();
        second_.print_size();
    }

private:
    SingleRead first_;
    SingleRead second_;
    size_t insert_size_;

};

inline std::ostream &operator<<(std::ostream &os, const PairedRead &read) {
    os << "Single read first=" << read.first() << " second=" << read.second() << std::endl;
    return os;
}

class PairedReadSeq {
public:
    typedef SingleReadSeq SingleReadT;
private:
    SingleReadSeq first_;
    SingleReadSeq second_;
    size_t insert_size_;

public:
    PairedReadSeq() : first_(), second_(), insert_size_(0) { }

    bool BinRead(std::istream &file, size_t is = 0) {
        first_.BinRead(file);
        second_.BinRead(file);

        insert_size_ = is - (size_t) first_.GetLeftOffset() - (size_t) second_.GetRightOffset();
        return !file.fail();
    }

    bool BinWrite(std::ostream &file, bool rc1 = false, bool rc2 = false) const {
        first_.BinWrite(file, rc1);
        second_.BinWrite(file, rc2);

        return !file.fail();
    }

    const SingleReadSeq &first() const {
        return first_;
    }

    const SingleReadSeq &second() const {
        return second_;
    }

    size_t insert_size() const {
        return insert_size_;
    }

    size_t distance() const {
        return insert_size_ - second_.size();
    }

    size_t gap() const {
        return insert_size_ - first_.size() - second_.size();
    }

    size_t size() const {
        return std::max(first_.size(), second_.size());
    }

    size_t nucl_count() const {
        return first_.size() + second_.size();
    }

    PairedReadSeq(const SingleReadSeq &first,
                  const SingleReadSeq &second,
                  size_t insert_size)
            : first_(first), second_(second), insert_size_(insert_size) { }

    const SingleReadSeq &operator[](size_t i) const {
        if (i == 0) {
            return first_;
        } else if (i == 1) {
            return second_;
        }
        VERIFY(false);
        return first_;
    }

    const PairedReadSeq operator!() const {
        return PairedReadSeq(!second_, !first_, insert_size_);
    }

};

inline std::ostream &operator<<(std::ostream &os, const PairedReadSeq &read) {
    os << "Paired read first=" << read.first() << " second=" << read.second() << std::endl;
    return os;
}

}
