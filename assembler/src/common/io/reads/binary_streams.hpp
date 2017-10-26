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

class BinaryFileSingleStream: public ReadStream<SingleReadSeq> {
    std::ifstream stream_;
    ReadStreamStat read_stat_;
    size_t current_;

    void Init() {
        stream_.clear();
        stream_.seekg(0);
        VERIFY(stream_.good());
        read_stat_.read(stream_);
        current_ = 0;
    }

public:

    BinaryFileSingleStream(const std::string& file_name_prefix, size_t file_num) {
        std::string fname;
        fname = file_name_prefix + "_" + std::to_string(file_num) + ".seq";
        stream_.open(fname.c_str(), std::ios_base::binary | std::ios_base::in);

        Init();
    }

    bool is_open() override {
        return stream_.is_open();
    }

    bool eof() override {
        return current_ == read_stat_.read_count;
    }

    BinaryFileSingleStream& operator>>(SingleReadSeq& read) override {
        read.BinRead(stream_);
        VERIFY(current_ < read_stat_.read_count);

        ++current_;
        return *this;
    }

    void close() override {
        current_ = 0;
        stream_.close();
    }

    void reset() override {
        Init();
    }

};

//returns FF oriented paired reads
class BinaryUnmergingPairedStream: public ReadStream<PairedReadSeq> {
    BinaryFileSingleStream stream_;
    size_t insert_size_;
    size_t read_length_;

    PairedReadSeq Convert(const SingleReadSeq &read) const {
        if (read.GetLeftOffset() >= read_length_ ||
                read.GetRightOffset() >= read_length_) {
            return PairedReadSeq();
        }
//        VERIFY(read_length_ >= read.GetLeftOffset() &&
//                       read_length_ >= read.GetRightOffset());

        const size_t left_length = std::min(read.size(),
                                      read_length_ - read.GetLeftOffset());
        const size_t right_length = std::min(read.size(),
                                      read_length_ - read.GetRightOffset());
        SingleReadSeq left(read.sequence().Subseq(0, left_length),
                           read.GetLeftOffset(), 0);
        SingleReadSeq right(read.sequence().Subseq(read.size() - right_length),
                            0, read.GetRightOffset());
        return PairedReadSeq(left, right, insert_size_);
    }

public:
    BinaryUnmergingPairedStream(const std::string& file_name_prefix,
                               size_t file_num,
                               size_t insert_size,
                               size_t read_length) :
            stream_(file_name_prefix, file_num),
            insert_size_(insert_size),
            read_length_(read_length) {
    }

    bool is_open() override {
        return stream_.is_open();
    }

    bool eof() override {
        return stream_.eof();
    }

    BinaryUnmergingPairedStream& operator>>(PairedReadSeq& read) override {
        SingleReadSeq single_read;
        stream_ >> single_read;
        read = Convert(single_read);
        return *this;
    }

    void close() override {
        stream_.close();
    }

    void reset() override {
        stream_.reset();
    }

};

class BinaryFilePairedStream: public ReadStream<PairedReadSeq> {
    std::ifstream stream_;
    size_t insert_size_;
    ReadStreamStat read_stat_;
    size_t current_;

    void Init() {
        stream_.clear();
        stream_.seekg(0);
        VERIFY(stream_.good());
        read_stat_.read(stream_);
        current_ = 0;
    }

public:

    BinaryFilePairedStream(const std::string& file_name_prefix, size_t file_num, size_t insert_szie): stream_(), insert_size_ (insert_szie) {
        std::string fname;
        fname = file_name_prefix + "_" + std::to_string(file_num) + ".seq";
        stream_.open(fname.c_str(), std::ios_base::binary | std::ios_base::in);

        Init();
    }

    bool is_open() override {
        return stream_.is_open();
    }

    bool eof() override {
        return current_ >= read_stat_.read_count;
    }

    BinaryFilePairedStream& operator>>(PairedReadSeq& read) override {
        read.BinRead(stream_, insert_size_);
        VERIFY(current_ < read_stat_.read_count);

        ++current_;
        return *this;
    }

    void close() override {
        current_ = 0;
        stream_.close();
    }


    void reset() override {
        Init();
    }

};

}
