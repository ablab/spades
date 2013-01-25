//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * binary_io.hpp
 *
 *  Created on: Apr 12, 2012
 *      Author: andrey
 */

#ifndef BINARY_IO_HPP_
#define BINARY_IO_HPP_

#include <fstream>

#include "verify.hpp"
#include "ireader.hpp"
#include "single_read.hpp"
#include "paired_read.hpp"

namespace io {

typedef io::IReader<io::SingleRead> SingleReadStream;
typedef io::IReader<io::PairedRead> PairedReadStream;


class BinaryWriter {

private:
    const std::string file_name_prefix_;

    size_t file_num_;

    std::vector<std::ofstream*> file_ds_;

    size_t buf_size_;


    template<class Read>
    void FlushBuffer(const std::vector<Read>& buffer, std::ostream& file, size_t from, size_t to) {
        for (size_t i = from; i < to; ++i) {
            buffer[i].BinWrite(file);
        }
    }

    template<class Read>
    void FlushBuffer(const std::vector<Read>& buffer, std::ostream& file) {
        FlushBuffer(buffer, file, 0, buffer.size());
    }

    template<class Read>
    ReadStat ToBinary(io::IReader<Read>& stream, size_t buf_size) {
        size_t buffer_reads = buf_size / (sizeof (Read) * 4);
        size_t reads_to_flush = buffer_reads * file_num_;

        std::vector< std::vector<Read> > buf(file_num_, std::vector<Read>(buffer_reads) );
        std::vector< ReadStat > read_stats(file_num_);
        std::vector< size_t > current_buf_sizes(file_num_, 0);
        size_t read_count = 0;

        for (size_t i = 0; i < file_num_; ++i) {
            file_ds_[i]->seekp(0);
            read_stats[i].write(*file_ds_[i]);
        }

        size_t buf_index;
        while (!stream.eof()) {
            buf_index = read_count % file_num_;

            stream >> buf[buf_index][current_buf_sizes[buf_index]];
            read_stats[buf_index].increase(buf[buf_index][current_buf_sizes[buf_index]]);

            ++current_buf_sizes[buf_index];
            VERBOSE_POWER(++read_count, " reads processed");

            if (read_count % reads_to_flush == 0) {
                for (size_t i = 0; i < file_num_; ++i) {
                    FlushBuffer(buf[i], *file_ds_[i]);
                    current_buf_sizes[i] = 0;
                }
            }
        }

        ReadStat result;
        for (size_t i = 0; i < file_num_; ++i) {
            buf[i].resize(current_buf_sizes[i]);
            FlushBuffer(buf[i], *file_ds_[i]);

            file_ds_[i]->seekp(0);
            read_stats[i].write(*file_ds_[i]);
            result.merge(read_stats[i]);
        }

        INFO(read_count << " reads written");
        return result;
    }


    template<class Read>
    ReadStat ToBinaryForThread(io::IReader<Read>& stream, size_t buf_size, size_t thread_num) {
        size_t buffer_reads = buf_size / (sizeof (Read) * 4);
        std::vector<Read> buf(buffer_reads);

        ReadStat stat;
        file_ds_[thread_num]->seekp(0);
        stat.write(*file_ds_[thread_num]);

        size_t current = 0;

        while (!stream.eof()) {
            stream >> buf[current];
            stat.increase(buf[current]);
            ++current;

            if (stat.read_count_ % buffer_reads == 0) {
                FlushBuffer(buf, *file_ds_[thread_num]);
                current = 0;
            }
        }

        buf.resize(current);
        FlushBuffer(buf, *file_ds_[thread_num]);

        file_ds_[thread_num]->seekp(0);
        stat.write(*file_ds_[thread_num]);

        return stat;
    }


public:

    BinaryWriter(const std::string& file_name_prefix, size_t file_num, size_t buf_size):
            file_name_prefix_(file_name_prefix), file_num_(file_num), file_ds_(), buf_size_(buf_size) {

        std::string fname;
        for (size_t i = 0; i < file_num_; ++i) {
            fname = file_name_prefix_ + "_" + ToString(i) + ".seq";
            file_ds_.push_back(new std::ofstream(fname, std::ios_base::binary));
        }
    }

    ~BinaryWriter() {
        for (size_t i = 0; i < file_num_; ++i) {
            if (file_ds_[i]->is_open()) {
                file_ds_[i]->close();
            }
            delete file_ds_[i];
        }
    }


    ReadStat ToBinary(io::IReader<io::SingleReadSeq>& stream) {
        return ToBinary(stream, buf_size_ / file_num_);
    }

    ReadStat ToBinary(io::IReader<io::SingleRead>& stream) {
        return ToBinary(stream, buf_size_ / file_num_);
    }

    ReadStat ToBinary(io::IReader<io::PairedReadSeq>& stream) {
        return ToBinary(stream, buf_size_ / (2 * file_num_));
    }

    ReadStat ToBinary(io::IReader<io::PairedRead>& stream) {
        return ToBinary(stream, buf_size_ / (2 * file_num_));
    }

    ReadStat ToBinaryForThread(io::IReader<io::SingleReadSeq>& stream, size_t thread_num) {
        return ToBinaryForThread(stream, buf_size_ / file_num_, thread_num);
    }

    ReadStat ToBinaryForThread(io::IReader<io::SingleRead>& stream, size_t thread_num) {
        return ToBinaryForThread(stream, buf_size_ / file_num_, thread_num);
    }

    ReadStat ToBinaryForThread(io::IReader<io::PairedReadSeq>& stream, size_t thread_num) {
        return ToBinaryForThread(stream, buf_size_ / (2 * file_num_), thread_num);
    }

    ReadStat ToBinaryForThread(io::IReader<io::PairedRead>& stream, size_t thread_num) {
        return ToBinaryForThread(stream, buf_size_ / (2 * file_num_), thread_num);
    }

//    template<class Read>
//    void WriteReads(std::vector<Read>& data) {
//        size_t chunk_size = data.size() / file_num_;
//        size_t last_chunk_size = chunk_size + data.size() % file_num_;
//
//        for (size_t i = 0; i < file_num_ - 1; ++i) {
//            file_ds_[i]->write((const char *) &chunk_size, sizeof(chunk_size));
//        }
//        file_ds_.back()->write((const char *) &last_chunk_size, sizeof(last_chunk_size));
//
//        size_t start_pos = 0;
//        for (size_t i = 0; i < file_num_ - 1; ++i, start_pos += chunk_size) {
//            FlushBuffer(data, *file_ds_[i], start_pos, start_pos + chunk_size);
//        }
//        FlushBuffer(data, file_ds_.back(), start_pos, data.size());
//    }
//
//    template<class Read>
//    void WriteSeparatedReads(std::vector< std::vector<Read> >& data) {
//        if (data.size() != file_num_) {
//            WARN("Cannot write reads, number of vectors is not equal to thread number");
//            return;
//        }
//
//        for (size_t i = 0; i < file_num_; ++i) {
//            size_t size = data[i].size();
//            file_ds_[i]->write((const char *) &size, sizeof(size));
//        }
//
//        for (size_t i = 0; i < file_num_; ++i) {
//            FlushBuffer(data[i], *file_ds_[i]);
//        }
//    }
};


template<class Read>
class PredictableIReader: public io::IReader<Read> {

public:
    virtual size_t size() const = 0;

};

// == Deprecated classes ==
// Use FileReadStream and InsertSizeModyfing instead

class SeqSingleReadStream: public io::PredictableIReader<io::SingleReadSeq> {

private:
    std::ifstream stream_;

    ReadStat read_stat_;

    size_t current_;

public:

    SeqSingleReadStream(const std::string& file_name_prefix, size_t file_num) {
        std::string fname;
        fname = file_name_prefix + "_" + ToString(file_num) + ".seq";
        stream_.open(fname.c_str(), std::ios_base::binary | std::ios_base::in);

        reset();
    }

    virtual ~SeqSingleReadStream() {
        if (stream_.is_open()) {
            stream_.close();
        }
    }

    virtual bool is_open() {
        return stream_.is_open();
    }

    virtual bool eof() {
        return current_ == read_stat_.read_count_;
    }

    virtual SeqSingleReadStream& operator>>(io::SingleReadSeq& read) {
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

    virtual ReadStat get_stat() const {
        return read_stat_;
    }

};


class SeqPairedReadStream: public io::PredictableIReader<io::PairedReadSeq> {

private:
    std::ifstream stream_;

    size_t insert_size_;

    ReadStat read_stat_;

    size_t current_;


public:

    SeqPairedReadStream(const std::string& file_name_prefix, size_t file_num, size_t insert_szie): stream_(), insert_size_ (insert_szie) {
        std::string fname;
        fname = file_name_prefix + "_" + ToString(file_num) + ".seq";
        stream_.open(fname.c_str(), std::ios_base::binary | std::ios_base::in);

        reset();
    }

    virtual ~SeqPairedReadStream() {
        if (stream_.is_open()) {
            stream_.close();
        }
    }

    virtual bool is_open() {
        return stream_.is_open();
    }

    virtual bool eof() {
        return current_ >= read_stat_.read_count_;
    }

    virtual SeqPairedReadStream& operator>>(io::PairedReadSeq& read) {
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

    ReadStat get_stat() const {
        ReadStat stat = read_stat_;
        stat.read_count_ *= 2;
        return stat;
    }
};


template <class Read>
class FileReadStream: public io::PredictableIReader<Read> {

private:
    std::ifstream stream_;

    ReadStat read_stat_;

    size_t current_;

public:

    FileReadStream(const std::string& file_name_prefix, size_t file_num) {
        std::string fname;
        fname = file_name_prefix + "_" + ToString(file_num) + ".seq";
        stream_.open(fname.c_str(), std::ios_base::binary | std::ios_base::in);

        reset();
    }

    virtual ~FileReadStream() {
        if (stream_.is_open()) {
            stream_.close();
        }
    }

    virtual bool is_open() {
        return stream_.is_open();
    }

    virtual bool eof() {
        return current_ == read_stat_.read_count_;
    }

    virtual SeqSingleReadStream& operator>>(Read& read) {
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

    virtual ReadStat get_stat() const {
        return read_stat_;
    }
};


template <class Read>
class ReadBufferedStream: public io::PredictableIReader<Read> {

private:
    std::vector<Read> * data_;

    ReadStat read_stat_;

    size_t current_;

public:

    ReadBufferedStream(io::PredictableIReader<Read>& stream) {
        read_stat_ = stream.get_stat();
        data_ = new std::vector<Read>(read_stat_.read_count_);

        size_t i = 0;
        while (!stream.eof()) {
            stream >> (*data_)[i++];
        }

        reset();
    }

    virtual ~ReadBufferedStream() {
        delete data_;
    }

    virtual bool is_open() {
        return true;
    }

    virtual bool eof() {
        return current_ == read_stat_.read_count_;
    }

    virtual ReadBufferedStream& operator>>(Read& read) {
        read = (*data_)[current_];
        VERIFY(current_ < read_stat_.read_count_);

        ++current_;
        return *this;
    }

    virtual void close() {
        current_ = 0;
    }

    virtual void reset() {
        current_ = 0;
    }

    virtual size_t size() const {
        return read_stat_.read_count_;
    }

    virtual ReadStat get_stat() const {
        return read_stat_;
    }
};


class SeqSingleReadStreamWrapper: public io::IReader<io::SingleReadSeq> {

private:
    io::IReader<io::PairedReadSeq>& stream_;

    PairedReadSeq current_read_;

    bool is_read_;

public:

    SeqSingleReadStreamWrapper(io::IReader<io::PairedReadSeq>& stream): stream_(stream), current_read_(), is_read_(false)  {
    }

    virtual ~SeqSingleReadStreamWrapper() {}

    virtual bool is_open() {
        return stream_.is_open();
    }

    virtual bool eof() {
        return stream_.eof() && !is_read_;
    }

    virtual SeqSingleReadStreamWrapper& operator>>(io::SingleReadSeq& read) {
        if (!is_read_) {
            stream_ >> current_read_;
            read = current_read_.first();
        } else {
            read = current_read_.second();
        }
        is_read_ = !is_read_;
        return *this;
    }

    virtual void close() {
        stream_.close();
    }

    virtual void reset() {
        stream_.reset();
        is_read_ = false;
    }

    virtual ReadStat get_stat() const {
        return stream_.get_stat();
    }
};



class CleanSeqSingleReadStreamWrapper: public io::IReader<io::SingleReadSeq> {

private:
    io::IReader<io::PairedReadSeq> * stream_;

    PairedReadSeq current_read_;

    bool is_read_;

public:

    CleanSeqSingleReadStreamWrapper(io::IReader<io::PairedReadSeq> * stream): stream_(stream), current_read_(), is_read_(false)  {
    }

    virtual ~CleanSeqSingleReadStreamWrapper() {
        delete stream_;
    }

    virtual bool is_open() {
        return stream_->is_open();
    }

    virtual bool eof() {
        return stream_->eof() && !is_read_;
    }

    virtual CleanSeqSingleReadStreamWrapper& operator>>(io::SingleReadSeq& read) {
        if (!is_read_) {
            stream_->operator >>(current_read_);
            read = current_read_.first();
        } else {
            read = current_read_.second();
        }
        is_read_ = !is_read_;
        return *this;
    }

    virtual void close() {
        stream_->close();
    }

    virtual void reset() {
        stream_->reset();
        is_read_ = false;
    }

    virtual ReadStat get_stat() const {
        return stream_->get_stat();
    }
};


class InsertSizeModifyingWrapper: public io::IReader<io::PairedReadSeq> {

private:
    io::IReader<io::PairedReadSeq>& stream_;

    size_t insert_size_;

public:

    InsertSizeModifyingWrapper(io::IReader<io::PairedReadSeq>& stream, size_t insert_szie): stream_(stream), insert_size_ (insert_szie) {
    }

    virtual ~InsertSizeModifyingWrapper() {
    }

    virtual bool is_open() {
        return stream_.is_open();
    }

    virtual bool eof() {
        return stream_.eof();
    }

    virtual InsertSizeModifyingWrapper& operator>>(io::PairedReadSeq& read) {
        stream_ >> read;
        read.inc_insert_size(insert_size_);
        return *this;
    }

    virtual void close() {
        stream_.close();
    }

    virtual void reset() {
        stream_.reset();
    }

    virtual ReadStat get_stat() const {
        return stream_.get_stat();
    }
};

}


#endif /* BINARY_IO_HPP_ */
