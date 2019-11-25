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
    typedef ReadStream<ReadType> ReadStreamT;

public:
    typedef ReadType ReadT;

    MultifileStream(ReadStreamList<ReadType> readers)
            : readers_{std::move(readers)}, current_reader_index_(0) {}

    MultifileStream(ReadStreamT reader_1) :
        current_reader_index_(0) {
        VERIFY(reader_1.is_open());
        readers_.push_back(std::move(reader_1));
    }

    MultifileStream(ReadStreamT reader_1, ReadStreamT reader_2) :
            current_reader_index_(0) {
        VERIFY(reader_1.is_open() && reader_2.is_open());
        readers_.push_back(std::move(reader_1));
        readers_.push_back(std::move(reader_2));
    }

    MultifileStream(MultifileStream<ReadType>&& multifile_stream) noexcept :
        readers_(std::move(multifile_stream.readers_)),
        current_reader_index_(std::exchange(multifile_stream.current_reader_index_, 0)) {}

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

    size_t size() const {
        return readers_.size();
    }

protected:
    ReadStreamList<ReadType> readers_;
    size_t current_reader_index_;
};


/**
 * ScopedMultifileStream is the stream that gets data from a number of streams,
 * given in a constructor and return ownership of streams to the initial state on destruction
 */
template<typename ReadType>
class ScopedMultifileStream : public MultifileStream<ReadType> {
    typedef ReadStream<ReadType> ReadStreamT;
public:
    typedef ReadType ReadT;

    ScopedMultifileStream(ReadStreamT& reader_1) : MultifileStream<ReadType>(std::move(reader_1)),
                                                   origin_stream_refs_({reader_1}) {}

    ScopedMultifileStream(ReadStreamT& reader_1, ReadStreamT& reader_2) :
        MultifileStream<ReadType>(std::move(reader_1), std::move(reader_2)),
        origin_stream_refs_({reader_1, reader_2}) {
    }

    ScopedMultifileStream(ScopedMultifileStream<ReadType>&& guard_multifile_stream) noexcept = default;

    ~ScopedMultifileStream() noexcept {
        for (size_t i = 0; i < origin_stream_refs_.size(); ++i) {
            origin_stream_refs_[i].get() = std::move(MultifileStream<ReadType>::readers_[i]);
        }
    }

private:
    std::vector<std::reference_wrapper<ReadStreamT>> origin_stream_refs_;
};


template<class ReadType>
ReadStream<ReadType> ScopedMultifileWrap(ReadStream <ReadType>& reader_1,
                                         ReadStream <ReadType>& reader_2) {
    return ScopedMultifileStream<ReadType>(reader_1, reader_2);
}

template<class ReadType>
ReadStream<ReadType> ScopedMultifileWrap(ReadStream <ReadType>& reader_1) {
    return ScopedMultifileStream<ReadType>(reader_1);
}

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
