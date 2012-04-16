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


class ReadsToBinaryConverter {

private:
    const std::string file_name_prefix_;

    size_t file_num_;

    std::vector<std::ofstream*> file_ds_;

    size_t buf_size_;

    template<class Read>
    void FlushBuffer(const std::vector<Read>& buffer, std::ostream& file) {
        for (auto iter = buffer.begin(); iter != buffer.end(); ++iter) {
            iter->BinWrite(file);
        }
    }

    template<class Read>
    void ToBinary(io::IReader<Read>& stream, size_t buf_size) {
        size_t read_count = 0;
        size_t reads_to_flush = buf_size * file_num_;

        std::vector< std::vector<Read> > buf(file_num_, std::vector<Read>(buf_size) );
        std::vector< size_t > buf_sizes(file_num_, 0);
        std::vector< size_t > current_buf_sizes(file_num_, 0);

        for (size_t i = 0; i < file_num_; ++i) {
            size_t read_num = 0;
            file_ds_[i]->write((const char *) &read_num, sizeof(read_num));
        }

        size_t buf_index;
        while (!stream.eof()) {
            buf_index = read_count % file_num_;

            stream >> buf[buf_index][current_buf_sizes[buf_index]];
            ++current_buf_sizes[buf_index];
            ++buf_sizes[buf_index];

            VERBOSE_POWER(++read_count, " reads processed");

            if (read_count == reads_to_flush) {
                for (size_t i = 0; i < file_num_; ++i) {
                    FlushBuffer(buf[i], *file_ds_[i]);
                    current_buf_sizes[i] = 0;
                }
            }
        }

        for (size_t i = 0; i < file_num_; ++i) {
            buf[i].resize(current_buf_sizes[i]);
            FlushBuffer(buf[i], *file_ds_[i]);

            file_ds_[i]->seekp(0);
            file_ds_[i]->write((const char *) &buf_sizes[i], sizeof(buf_sizes[i]));
        }

        INFO(read_count << " reads converted");
    }

public:

    ReadsToBinaryConverter(const std::string& file_name_prefix, size_t file_num, size_t buf_size):
            file_name_prefix_(file_name_prefix), file_num_(file_num), file_ds_(), buf_size_(buf_size) {
        std::string fname;
        for (size_t i = 0; i < file_num_; ++i) {
            fname = file_name_prefix_ + "_" + ToString(i) + ".seq";
            file_ds_.push_back(new std::ofstream());
            file_ds_.back()->open(fname.c_str(), std::ios_base::binary | std::ios_base::out);
        }
    }

    ~ReadsToBinaryConverter() {
        for (size_t i = 0; i < file_num_; ++i) {
            if (file_ds_[i]->is_open()) {
                file_ds_[i]->close();
            }
            delete file_ds_[i];
        }
    }

    void ToBinary(SingleReadStream& stream) {
        ToBinary(stream, buf_size_ / file_num_);
    }

    void ToBinary(PairedReadStream& stream) {
        ToBinary(stream, (buf_size_ / file_num_ * 2));
    }

};


class SeqSingleReadStream: public io::IReader<io::SingleReadSeq> {

private:
    std::ifstream stream_;

    size_t read_num_;

    size_t current_;

public:

    SeqSingleReadStream(const std::string& file_name_prefix, size_t file_num) {
        std::string fname;
        fname = file_name_prefix + "_" + ToString(file_num) + ".seq";
        stream_.open(fname.c_str(), std::ios_base::binary | std::ios_base::in);

        stream_.read((char *) &read_num_, sizeof(read_num_));
        current_ = 0;
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
        return current_ >= read_num_;
    }

    virtual SeqSingleReadStream& operator>>(io::SingleReadSeq& read) {
        read.BinRead(stream_);
        ++current_;
        return *this;
    }

    virtual void close() {
        current_ = 0;
        stream_.close();
    }

    virtual void reset() {
        current_ = 0;
        stream_.seekg(0);
    }
};


class SeqPairedReadStream: public io::IReader<io::PairedReadSeq> {

private:
    std::ifstream stream_;

    size_t insert_size_;

    size_t read_num_;

    size_t current_;


public:

    SeqPairedReadStream(const std::string& file_name_prefix, size_t file_num, size_t insert_szie): stream_(), insert_size_ (insert_szie) {
        std::string fname;
        fname = file_name_prefix + "_" + ToString(file_num) + ".seq";
        stream_.open(fname.c_str(), std::ios_base::binary | std::ios_base::in);

        stream_.read((char *) &read_num_, sizeof(read_num_));
        current_ = 0;
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
        return current_ >= read_num_;
    }

    virtual SeqPairedReadStream& operator>>(io::PairedReadSeq& read) {
        read.BinRead(stream_, insert_size_);
        ++current_;
        return *this;
    }

    virtual void close() {
        current_ = 0;
        stream_.close();
    }

    virtual void reset() {
        current_ = 0;
        stream_.seekg(0);
    }
};

class SeqSingleReadStreamWrapper: public io::IReader<io::SingleReadSeq> {

private:
    SeqPairedReadStream& stream_;

    PairedReadSeq current_read_;

    bool is_read_;

public:

    SeqSingleReadStreamWrapper(SeqPairedReadStream& stream): stream_(stream), current_read_(), is_read_(false)  {
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

};

}


#endif /* BINARY_IO_HPP_ */
