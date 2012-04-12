/*
 * binary_io.hpp
 *
 *  Created on: Apr 12, 2012
 *      Author: andrey
 */

#ifndef BINARY_IO_HPP_
#define BINARY_IO_HPP_

#include <fstream>

#include <verify.hpp>
#include <ireader.hpp>
#include <single_read.hpp>
#include <paired_read.hpp>

namespace io {

typedef io::IReader<io::SingleRead> SingleReadStream;
typedef io::IReader<io::PairedRead> PairedReadStream;

class ReadsToBinaryConverter {

private:
    const std::string file_name_prefix_;
    size_t file_num_;
    std::vector<std::ofstream> file_ds_;

public:

    ReadsToBinaryConverter(const std::string& file_name_prefix, size_t file_num): file_name_prefix_(file_name_prefix), file_num_(file_num), file_ds_(file_num, std::ofstream()) {
        std::string fname;
        for (size_t i = 0; i < file_num_; ++i) {
            fname = file_name_prefix_ + "_" + ToString(i) + ".seq";
            file_ds_[i].open(fname.c_str(), std::ios_base::binary | std::ios_base::out);
        }
    }

    ~ReadsToBinaryConverter() {
        for (size_t i = 0; i < file_num_; ++i) {
            if (file_ds_[i].is_open())
                file_ds_[i].close();
        }
    }

    void ToBinary(SingleReadStream& stream) {
        io::SingleRead sr;
        size_t read_count = 0;

        while (!stream.eof()) {
            stream >> sr;
            sr.sequence().BinWrite(file_ds_[read_count % file_num_]);
            ++read_count;
        }
    }

    void ToBinary(PairedReadStream& stream) {
        io::PairedRead sr;
        size_t read_count = 0;
        int index = 0;

        while (!stream.eof()) {
            index = read_count % file_num_;

            stream >> sr;
            sr.first().sequence().BinWrite(file_ds_[index]);
            sr.second().sequence().BinWrite(file_ds_[index]);
            ++read_count;
        }
    }
};


class SeqSingleReadStream: public IReader<SingleReadSeq> {

private:
    std::ifstream& stream_;

public:

    SeqSingleReadStream(std::ifstream& stream): stream_(stream) {
    }

    virtual ~SeqSingleReadStream() {}

    virtual bool is_open() {
        return stream_.is_open();
    }

    virtual bool eof() {
        return stream_.eof();
    }

    virtual SeqSingleReadStream& operator>>(SingleReadSeq& read) {
        return SingleReadSeq(stream_);
    }

    virtual void close() {
        stream_.close();
    }

    virtual void reset() = 0;
};


class ReadsFromBinaryConverter {

private:
    const std::string file_name_prefix_;
    size_t file_num_;
    std::vector<std::ifstream> file_ds_;

public:

    ReadsFromBinaryConverter(const std::string& file_name_prefix, size_t file_num): file_name_prefix_(file_name_prefix), file_num_(file_num), file_ds_(file_num, std::ifstream()) {
        std::string fname;
        for (size_t i = 0; i < file_num_; ++i) {
            fname = file_name_prefix_ + "_" + ToString(i) + ".seq";
            file_ds_[i].open(fname.c_str(), std::ios_base::binary | std::ios_base::in);
        }
    }

    ~ReadsFromBinaryConverter() {
        for (size_t i = 0; i < file_num_; ++i) {
            if (file_ds_[i].is_open())
                file_ds_[i].close();
        }
    }

    std::ifstream& getRawStream(int index) const {
        VERIFY(index < file_num_);
        return file_ds_[index];
    }

};

}


#endif /* BINARY_IO_HPP_ */
