//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "utils/verify.hpp"
#include "ireader.hpp"
#include "single_read.hpp"
#include "paired_read.hpp"
#include "binary_converter.hpp"

#include <fstream>

namespace io {

template<typename SeqT>
class BinaryFileStream: public ReadStream<SeqT> {

protected:
    std::ifstream stream_;

    virtual bool ReadImpl(SeqT &read) = 0;

private:
    size_t offset_, count_, current_;

    void Init() {
        stream_.clear();
        stream_.seekg(offset_);
        VERIFY(stream_.good());
        current_ = 0;
    }

public:
    BinaryFileStream(const std::string &file_name_prefix, size_t total_num, size_t portion_num)
            : offset_(sizeof(ReadStreamStat)) {
        INFO("Preparing binary stream #" << portion_num << "/" << total_num);
        VERIFY(portion_num < total_num);
        std::string fname = file_name_prefix + ".seq";
        stream_.open(fname, std::ios_base::binary | std::ios_base::in);
        ReadStreamStat stat;
        stat.read(stream_);

        //Determining the number of N-read chunks and actual number of reads
        std::string offset_name = file_name_prefix + ".off";
        const size_t chunk_count = fs::filesize(offset_name) / sizeof(size_t);
        const size_t chunk_num = chunk_count * portion_num / total_num;

        if (chunk_count) {
            //Calculating the absolute offset in the reads file
            std::ifstream offset_stream(offset_name, std::ios_base::binary | std::ios_base::in);
            offset_stream.seekg(chunk_num * sizeof(size_t));
            offset_stream.read(reinterpret_cast<char *>(&offset_), sizeof(offset_));
        }

        //Calculating the size of our portion (be careful at the file end)
        const size_t start_num = chunk_num * BinaryWriter::CHUNK;
        const size_t start_next = chunk_count * (portion_num + 1) / total_num * BinaryWriter::CHUNK;
        count_ = std::min(stat.read_count, start_next) - start_num;

        INFO("Reads " << start_num << "-" << start_num + count_ << "/" << stat.read_count
                       << " at 0x" << std::ios_base::hex << offset_);

        Init();
    }

    BinaryFileStream<SeqT>& operator>>(SeqT &read) override {
        VERIFY(current_++ < count_);
        VERIFY(ReadImpl(read));
        return *this;
    }

    bool is_open() override {
        return stream_.is_open();
    }

    bool eof() override {
        return current_ >= count_;
    }

    void close() override {
        current_ = 0;
        stream_.close();
    }

    void reset() override {
        Init();
    }

};

class BinaryFileSingleStream : public BinaryFileStream<SingleReadSeq>  {

protected:
    bool ReadImpl(SingleReadSeq &read) override {
        return read.BinRead(stream_);
    }

public:
    BinaryFileSingleStream(const std::string &file_name_prefix, size_t total_num, size_t portion_num)
            : BinaryFileStream(file_name_prefix, total_num, portion_num) {
    }
};

class BinaryFilePairedStream: public BinaryFileStream<PairedReadSeq> {
    size_t insert_size_;

protected:
    bool ReadImpl(PairedReadSeq& read) override {
        return read.BinRead(stream_, insert_size_);
    }

public:
    BinaryFilePairedStream(const std::string &file_name_prefix, size_t insert_size,
                           size_t total_num, size_t portion_num)
            : BinaryFileStream(file_name_prefix, total_num, portion_num), insert_size_ (insert_size) {
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

        const size_t left_length = std::min(read.size(), read_length_ - read.GetLeftOffset());
        const size_t right_length = std::min(read.size(), read_length_ - read.GetRightOffset());
        SingleReadSeq left(read.sequence().Subseq(0, left_length), read.GetLeftOffset(), 0);
        SingleReadSeq right(read.sequence().Subseq(read.size() - right_length), 0, read.GetRightOffset());
        return PairedReadSeq(left, right, insert_size_);
    }

public:
    BinaryUnmergingPairedStream(const std::string& file_name_prefix, size_t insert_size, size_t read_length,
                                size_t total_num, size_t portion_num) :
            stream_(file_name_prefix, total_num, portion_num),
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

}
