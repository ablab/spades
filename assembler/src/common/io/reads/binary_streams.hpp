//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <fstream>

#include "utils/verify.hpp"
#include "ireader.hpp"
#include "single_read.hpp"
#include "paired_read.hpp"

namespace io {

// == Deprecated classes ==
// Use FileReadStream and InsertSizeModyfing instead

class BinaryFileSingleStream: public PredictableReadStream<SingleReadSeq> {
private:
    std::ifstream stream_;
    ReadStreamStat read_stat_;
    size_t current_;

public:

    BinaryFileSingleStream(const std::string& file_name_prefix, size_t file_num) {
        std::string fname;
        fname = file_name_prefix + "_" + ToString(file_num) + ".seq";
        stream_.open(fname.c_str(), std::ios_base::binary | std::ios_base::in);

        reset();
    }

    virtual bool is_open() {
        return stream_.is_open();
    }

    virtual bool eof() {
        return current_ == read_stat_.read_count_;
    }

    virtual BinaryFileSingleStream& operator>>(SingleReadSeq& read) {
        read.BinRead(stream_);
        VERIFY(current_ < read_stat_.read_count_);

        ++current_;
        return *this;
    }

    virtual void close() {
        current_ = 0;
        stream_.close();
    }

    virtual void reset() {
        stream_.clear();
        stream_.seekg(0);
        VERIFY(stream_.good());
        read_stat_.read(stream_);
        current_ = 0;
    }

    virtual size_t size() const {
        return read_stat_.read_count_;
    }

    virtual ReadStreamStat get_stat() const {
        return read_stat_;
    }

};

class BinaryFilePairedStream: public PredictableReadStream<PairedReadSeq> {

private:
    std::ifstream stream_;

    size_t insert_size_;

    ReadStreamStat read_stat_;

    size_t current_;


public:

    BinaryFilePairedStream(const std::string& file_name_prefix, size_t file_num, size_t insert_szie): stream_(), insert_size_ (insert_szie) {
        std::string fname;
        fname = file_name_prefix + "_" + ToString(file_num) + ".seq";
        stream_.open(fname.c_str(), std::ios_base::binary | std::ios_base::in);

        reset();
    }

    virtual bool is_open() {
        return stream_.is_open();
    }

    virtual bool eof() {
        return current_ >= read_stat_.read_count_;
    }

    virtual BinaryFilePairedStream& operator>>(PairedReadSeq& read) {
        read.BinRead(stream_, insert_size_);
        VERIFY(current_ < read_stat_.read_count_);

        ++current_;
        return *this;
    }

    virtual void close() {
        current_ = 0;
        stream_.close();
    }


    virtual void reset() {
        stream_.clear();
        stream_.seekg(0);
        VERIFY(stream_.good());
        read_stat_.read(stream_);
        current_ = 0;
    }

    virtual size_t size() const {
        return read_stat_.read_count_;
    }

    ReadStreamStat get_stat() const {
        ReadStreamStat stat = read_stat_;
        stat.read_count_ *= 2;
        return stat;
    }
};

}
