//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "read_stream_vector.hpp"
#include <vector>

namespace io {

/**
 * MultifileReader is the stream that gets data from number of files,
 * given in a constructor.
 */
template<typename ReadType>
class MultifileStream: public ReadStream<ReadType> {
    typedef ReadStream<ReadType> StreamT;
    typedef std::shared_ptr<StreamT> ReadStreamPtrT;
public:
    MultifileStream(const ReadStreamList<ReadType>& readers) :
        readers_(readers), current_reader_index_(0) {
    }

    MultifileStream(ReadStreamPtrT reader_1, ReadStreamPtrT reader_2) :
            current_reader_index_(0) {
        VERIFY(reader_1->is_open() && reader_2->is_open());
        readers_.push_back(reader_1);
        readers_.push_back(reader_2);
    }

    /* virtual */
    bool is_open() {
        return (readers_.size() > 0) && readers_[0].is_open();
    }

    /* virtual */
    bool eof() {
        while ((current_reader_index_ < readers_.size()) && readers_[current_reader_index_].eof()) {
            ++current_reader_index_;
        }
        return current_reader_index_ == readers_.size();
    }

    /* virtual */
    MultifileStream& operator>>(ReadType& read) {
        if (!eof()) {
            readers_[current_reader_index_] >> read;
        }
        return (*this);
    }

    /* virtual */
    void close() {
        readers_.close();
    }

    /* virtual */
    void reset() {
        readers_.reset();
        current_reader_index_ = 0;
    }

    /* virtual */
    ReadStreamStat get_stat() const {
        return readers_.get_stat();
    }

private:
    ReadStreamList<ReadType> readers_;
    size_t current_reader_index_;
};

template<class ReadType>
std::shared_ptr<ReadStream<ReadType>> MultifileWrap(std::shared_ptr<ReadStream<ReadType>> reader_1,
                                                  std::shared_ptr<ReadStream<ReadType>> reader_2) {
    return std::make_shared<MultifileStream<ReadType>>(reader_1, reader_2);
}

template<class ReadType>
std::shared_ptr<ReadStream<ReadType>> MultifileWrap(const ReadStreamList<ReadType>& readers) {
    return std::make_shared<MultifileStream<ReadType>>(readers);
}

template<class ReadType>
ReadStreamList<ReadType> WrapPairsInMultifiles(ReadStreamList<ReadType> readers_1,
                                           ReadStreamList<ReadType> readers_2) {
    VERIFY(readers_1.size() == readers_2.size());
    ReadStreamList<ReadType> answer;
    for (size_t i = 0; i < readers_1.size(); ++i) {
        answer.push_back(MultifileWrap<ReadType>(readers_1.ptr_at(i), readers_2.ptr_at(i)));
    }
    return answer;
}

}
