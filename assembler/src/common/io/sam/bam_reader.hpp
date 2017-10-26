//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************
//todo rename to reader
#pragma once

#include "io/reads/ireader.hpp"
#include "io/reads/single_read.hpp"

#include <bamtools/api/BamReader.h>

namespace io {
class BamRead : public BamTools::BamAlignment {
public:
    BamRead() { }

    BamRead(const BamTools::BamAlignment &other)
            : BamTools::BamAlignment(other) { }

    const std::string &name() const {
        return Name;
    }

    size_t size() const {
        return Length;
    }

    size_t nucl_count() const {
        return size();
    }

    const std::string &GetSequenceString() const {
        return QueryBases;
    }

    std::string GetPhredQualityString() const {
        return Qualities;
    }

    operator io::SingleRead() {
        // not including quality is intentional:
        // during read correction bases might be inserted/deleted,
        // and base qualities for them are not calculated
        return io::SingleRead(name(), GetSequenceString());
    }

    char operator[](size_t i) const {
        VERIFY(is_nucl(QueryBases[i]));
        return dignucl(QueryBases[i]);
    }
};

class UnmappedBamStream : public ReadStream<BamRead> {
public:
    UnmappedBamStream(const std::string &filename)
            : filename_(filename) {
        open();
    }

    virtual ~UnmappedBamStream() { }

    bool is_open() { return is_open_; }

    bool eof() { return eof_; }

    UnmappedBamStream &operator>>(BamRead &read) {
        if (!is_open_ || eof_)
            return *this;

        read = seq_;
        eof_ = (false == reader_.GetNextAlignment(seq_));

        return *this;
    }

    void close() {
        reader_.Close();
        is_open_ = false;
        eof_ = true;
    }

    void reset() {
        close();
        open();
    }

private:
    BamTools::BamReader reader_;
    BamTools::BamAlignment seq_;
    std::string filename_;
    bool is_open_;
    bool eof_;

    void open() {
        reader_.Open(filename_);
        is_open_ = true;

        eof_ = (false == reader_.GetNextAlignment(seq_));
    }

};
}
