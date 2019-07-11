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
class MultifileStream {
    typedef ReadStream<ReadType> ReadStreamPtrT;

public:
    MultifileStream(ReadStreamList<ReadType> readers)
            : readers_{std::move(readers)}, current_reader_index_(0) {}

    MultifileStream(ReadStreamPtrT reader_1, ReadStreamPtrT reader_2) :
            current_reader_index_(0) {
        VERIFY(reader_1.is_open() && reader_2.is_open());
        readers_.push_back(std::move(reader_1));
        readers_.push_back(std::move(reader_2));
    }

    bool is_open() {
        return (readers_.size() > 0) && readers_[0].is_open();
    }

    bool eof() {
        while ((current_reader_index_ < readers_.size()) && readers_[current_reader_index_].eof()) {
            ++current_reader_index_;
        }
        return current_reader_index_ == readers_.size();
    }

    MultifileStream& operator>>(ReadType& read) {
        if (!eof()) {
            readers_[current_reader_index_] >> read;
        }
        return (*this);
    }

    void close() {
        readers_.close();
    }

    void reset() {
        readers_.reset();
        current_reader_index_ = 0;
    }

private:
    ReadStreamList<ReadType> readers_;
    size_t current_reader_index_;
};

template<class ReadType>
ReadStream<ReadType> MultifileWrap(ReadStream<ReadType> reader_1,
                                      ReadStream<ReadType> reader_2) {
    return MultifileStream<ReadType>(std::move(reader_1), std::move(reader_2));
}

template<class ReadType>
ReadStream<ReadType> MultifileWrap(ReadStreamList<ReadType> readers) {
    return MultifileStream<ReadType>(std::move(readers));
}

template<class ReadType>
ReadStreamList<ReadType> WrapPairsInMultifiles(ReadStreamList<ReadType> readers_1,
                                               ReadStreamList<ReadType> readers_2) {
    VERIFY(readers_1.size() == readers_2.size());
    ReadStreamList<ReadType> answer;
    for (size_t i = 0; i < readers_1.size(); ++i) {
        answer.push_back(MultifileWrap<ReadType>(std::move(readers_1[i]), std::move(readers_2[i])));
    }
    return answer;
}

}
