/*
 * binary_streams.hpp
 *
 *  Created on: Mar 28, 2013
 *      Author: andrey
 */

#ifndef BINARY_STREAMS_HPP_
#define BINARY_STREAMS_HPP_

#include <fstream>

#include "verify.hpp"
#include "ireader.hpp"
#include "single_read.hpp"
#include "paired_read.hpp"

namespace io {

typedef io::IReader<io::SingleRead> SingleReadStream;
typedef io::IReader<io::PairedRead> PairedReadStream;



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


typedef io::IReader<io::SingleReadSeq> SequenceSingleReadStream;
typedef io::IReader<io::PairedReadSeq> SequencePairedReadStream;

inline std::vector< SequenceSingleReadStream* > apply_single_wrappers(bool followed_by_rc,
        std::vector< SequenceSingleReadStream* >& single_readers,
        std::vector< SequencePairedReadStream* > * paired_readers = 0) {

    VERIFY(single_readers.size() != 0);
    size_t size = single_readers.size();
    std::vector<SequenceSingleReadStream*> raw_readers(size);

    if (paired_readers != 0) {
        VERIFY(single_readers.size() == paired_readers->size());

        for (size_t i = 0; i < size; ++i) {
            SequenceSingleReadStream * single_stream = single_readers.at(i);
            SequencePairedReadStream * paired_stream = paired_readers->at(i);
            io::CleanSeqSingleReadStreamWrapper * single_wrapper = new io::CleanSeqSingleReadStreamWrapper(paired_stream);

            raw_readers[i] = new io::MultifileReader<io::SingleReadSeq>(*single_wrapper, *single_stream, true);
        }
    }
    else {
       for (size_t i = 0; i < size; ++i) {
           raw_readers[i] = single_readers.at(i);
       }
    }

    if (followed_by_rc) {
        std::vector<SequenceSingleReadStream*> rc_readers(size);
        for (size_t i = 0; i < size; ++i) {
            rc_readers[i] = new io::CleanRCReaderWrapper<io::SingleReadSeq>(raw_readers[i]);
        }
        return rc_readers;
    } else {
        return raw_readers;
    }
}


inline std::vector< SequencePairedReadStream* > apply_paired_wrappers(bool followed_by_rc,
        std::vector< SequencePairedReadStream* >& paired_readers) {

    VERIFY(paired_readers.size() != 0);
    size_t size = paired_readers.size();

    if (followed_by_rc) {
        std::vector<SequencePairedReadStream*> rc_readers(size);
        for (size_t i = 0; i < size; ++i) {
            rc_readers[i] = new io::CleanRCReaderWrapper<io::PairedReadSeq>(paired_readers[i]);
        }
        return rc_readers;
    } else {
        return paired_readers;
    }
}



}


#endif /* BINARY_STREAMS_HPP_ */
