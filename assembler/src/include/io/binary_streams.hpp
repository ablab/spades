#pragma once

#include <fstream>

#include "verify.hpp"
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


//template <class Read>
//class FileReadStream: public io::PredictableIReader<Read> {
//
//private:
//    std::ifstream stream_;
//
//    ReadStat read_stat_;
//
//    size_t current_;
//
//public:
//
//    FileReadStream(const std::string& file_name_prefix, size_t file_num) {
//        std::string fname;
//        fname = file_name_prefix + "_" + ToString(file_num) + ".seq";
//        stream_.open(fname.c_str(), std::ios_base::binary | std::ios_base::in);
//
//        reset();
//    }
//
//    virtual ~FileReadStream() {
//        if (stream_.is_open()) {
//            stream_.close();
//        }
//    }
//
//    virtual bool is_open() {
//        return stream_.is_open();
//    }
//
//    virtual bool eof() {
//        return current_ == read_stat_.read_count_;
//    }
//
//    virtual FileReadStream& operator>>(Read& read) {
//        read.BinRead(stream_);
//        VERIFY(current_ < read_stat_.read_count_);
//
//        ++current_;
//        return *this;
//    }
//
//    virtual void close() {
//        current_ = 0;
//        stream_.close();
//    }
//
//    virtual void reset() {
//        stream_.clear();
//        stream_.seekg(0);
//        VERIFY(stream_.good());
//        read_stat_.read(stream_);
//        current_ = 0;
//    }
//
//    virtual size_t size() const {
//        return read_stat_.read_count_;
//    }
//
//    virtual ReadStat get_stat() const {
//        return read_stat_;
//    }
//};

//template <class Read>
//class ReadBufferedStream: public io::PredictableIReader<Read> {
//
//private:
//    std::vector<Read> * data_;
//
//    ReadStat read_stat_;
//
//    size_t current_;
//
//public:
//
//    ReadBufferedStream(io::PredictableIReader<Read>& stream) {
//        read_stat_ = stream.get_stat();
//        data_ = new std::vector<Read>(read_stat_.read_count_);
//
//        size_t i = 0;
//        while (!stream.eof()) {
//            stream >> (*data_)[i++];
//        }
//
//        reset();
//    }
//
//    virtual ~ReadBufferedStream() {
//        delete data_;
//    }
//
//    virtual bool is_open() {
//        return true;
//    }
//
//    virtual bool eof() {
//        return current_ == read_stat_.read_count_;
//    }
//
//    virtual ReadBufferedStream& operator>>(Read& read) {
//        read = (*data_)[current_];
//        VERIFY(current_ < read_stat_.read_count_);
//
//        ++current_;
//        return *this;
//    }
//
//    virtual void close() {
//        current_ = 0;
//    }
//
//    virtual void reset() {
//        current_ = 0;
//    }
//
//    virtual size_t size() const {
//        return read_stat_.read_count_;
//    }
//
//    virtual ReadStat get_stat() const {
//        return read_stat_;
//    }
//};

//class SeqSingleReadStreamWrapper: public Reader<SingleReadSeq> {
//
//private:
//    io::IReader<io::PairedReadSeq>& stream_;
//
//    PairedReadSeq current_read_;
//
//    bool is_read_;
//
//public:
//
//    SeqSingleReadStreamWrapper(io::IReader<io::PairedReadSeq>& stream): stream_(stream), current_read_(), is_read_(false)  {
//    }
//
//    virtual ~SeqSingleReadStreamWrapper() {}
//
//    virtual bool is_open() {
//        return stream_.is_open();
//    }
//
//    virtual bool eof() {
//        return stream_.eof() && !is_read_;
//    }
//
//    virtual SeqSingleReadStreamWrapper& operator>>(io::SingleReadSeq& read) {
//        if (!is_read_) {
//            stream_ >> current_read_;
//            read = current_read_.first();
//        } else {
//            read = current_read_.second();
//        }
//        is_read_ = !is_read_;
//        return *this;
//    }
//
//    virtual void close() {
//        stream_.close();
//    }
//
//    virtual void reset() {
//        stream_.reset();
//        is_read_ = false;
//    }
//
//    virtual ReadStat get_stat() const {
//        return stream_.get_stat();
//    }
//};

//class InsertSizeModifyingWrapper: public io::IReader<io::PairedReadSeq> {
//
//private:
//    io::IReader<io::PairedReadSeq>& stream_;
//
//    size_t insert_size_;
//
//public:
//
//    InsertSizeModifyingWrapper(io::IReader<io::PairedReadSeq>& stream, size_t insert_szie): stream_(stream), insert_size_ (insert_szie) {
//    }
//
//    virtual ~InsertSizeModifyingWrapper() {
//    }
//
//    virtual bool is_open() {
//        return stream_.is_open();
//    }
//
//    virtual bool eof() {
//        return stream_.eof();
//    }
//
//    virtual InsertSizeModifyingWrapper& operator>>(io::PairedReadSeq& read) {
//        stream_ >> read;
//        read.inc_insert_size(insert_size_);
//        return *this;
//    }
//
//    virtual void close() {
//        stream_.close();
//    }
//
//    virtual void reset() {
//        stream_.reset();
//    }
//
//    virtual ReadStat get_stat() const {
//        return stream_.get_stat();
//    }
//};

}
