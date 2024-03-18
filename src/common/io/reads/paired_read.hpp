//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "sequence/sequence.hpp"
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

    bool BinWrite(std::ostream &file, bool rc1 = false, bool rc2 = false,
                  uint64_t tag1 = 0, uint64_t tag2 = 0) const {
        first_.BinWrite(file, rc1, tag1);
        second_.BinWrite(file, rc2, tag2);

        return !file.fail();
    }

    UniversalPairedRead() : first_(), second_(), insert_size_(0) { }

    UniversalPairedRead(SingleReadT first,
                        SingleReadT second,
                        size_t insert_size)
            : first_{std::move(first)}, second_{std::move(second)},
              insert_size_(insert_size) { }


private:
    SingleReadT first_;
    SingleReadT second_;
    size_t insert_size_;
};

typedef UniversalPairedRead<SingleRead> PairedRead;

class TellSeqRead : public PairedRead {
  public:
    const auto &aux() const { return aux_; }
    auto &aux() { return aux_; }

    const TellSeqRead operator!() const {
        return TellSeqRead(!second(), !first(), orig_insert_size(), aux_);
    }


    TellSeqRead() = default;

    TellSeqRead(SingleReadT first, SingleReadT second,
                size_t insert_size,
                SingleReadT aux)
            : PairedRead(std::move(first), std::move(second), insert_size),
              aux_(aux) {}

  private:
    SingleReadT aux_;
};

inline std::ostream &operator<<(std::ostream &os, const PairedRead &read) {
    os << "Single read first=" << read.first() << " second=" << read.second() << std::endl;
    return os;
}

inline std::ostream &operator<<(std::ostream &os, const TellSeqRead &read) {
    os <<  "Single read first=" << read.first() << " second=" << read.second() << " aux=" << read.aux() << std::endl;
    return os;
}

typedef UniversalPairedRead<SingleReadSeq> PairedReadSeq;

inline std::ostream &operator<<(std::ostream &os, const PairedReadSeq &read) {
    os << "Paired read first=" << read.first() << " second=" << read.second() << std::endl;
    return os;
}

}
