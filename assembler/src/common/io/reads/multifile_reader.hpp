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

    void get_stream(size_t i, ReadStreamT& reader) {
        VERIFY(i < readers_.size());
        reader = std::move(readers_[i]);
    }

    void split_streams(ReadStreamT& reader_1, ReadStreamT& reader_2) {
        VERIFY(readers_.size() == 2);
        reader_1 = std::move(readers_[0]);
        reader_2 = std::move(readers_[1]);
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

template<class ReadType>
class MultifileReadStreamList {
public:
    typedef ReadType ReadT;
    typedef MultifileStream<ReadType> ReaderT;
    typedef std::vector<ReaderT> ReadersT;
    using iterator = typename ReadersT::iterator;
    using const_iterator = typename ReadersT::const_iterator;

private:
    ReadersT streams_;

public:
    MultifileReadStreamList(ReadStreamList<ReadType> streams1,
                            ReadStreamList<ReadType> streams2) {
        if (streams2.size() == 0) {
            for (size_t i = 0; i < streams1.size(); ++i) {
                streams_.push_back(MultifileStream<ReadType>(std::move(streams1[i])));
            }
            return;
        }

        VERIFY(streams1.size() == streams2.size());
        for (size_t i = 0; i < streams1.size(); ++i) {
            streams_.push_back(MultifileStream<ReadType>(std::move(streams1[i]), std::move(streams2[i])));
        }
    }

    MultifileReadStreamList() = default;
    MultifileReadStreamList(const MultifileReadStreamList&) = delete;
    MultifileReadStreamList(MultifileReadStreamList&&) = default;
    MultifileReadStreamList &operator=(const MultifileReadStreamList&) = delete;
    MultifileReadStreamList &operator=(MultifileReadStreamList&&) = default;

    void split_streams(ReadStreamList<ReadType>& reader_1, ReadStreamList<ReadType>& reader_2) {
        if (streams_.size() > 0 && streams_[0].size() == 1) {
            for (size_t i = 0; i < streams_.size(); ++i) {
                if (reader_1.size() <= i) {
                    reader_1.push_back(ReadStream<ReadType>());
                }

                streams_[i].get_stream(0, reader_1[i]);
            }
            return;
        }

        for (size_t i = 0; i < streams_.size(); ++i) {
            if (reader_1.size() <= i) {
                reader_1.push_back(ReadStream<ReadType>());
            }
            if (reader_2.size() <= i) {
                reader_2.push_back(ReadStream<ReadType>());
            }

            streams_[i].split_streams(reader_1[i], reader_2[i]);
        }
    }

    ReaderT &operator[](size_t i) { return streams_[i]; }
    const ReaderT &operator[](size_t i) const { return streams_[i]; }

    ReaderT &back() {
        return streams_.back();
    }

    size_t size() const {
        return streams_.size();
    }

    bool eof()  {
        for (size_t i = 0; i < streams_.size(); ++i) {
            if (!streams_[i].eof()) {
                return false;
            }
        }
        return true;
    }

    iterator begin() { return streams_.begin(); }
    iterator end() { return iterator(streams_.end()); }
    const_iterator begin() const { return streams_.begin(); }
    const_iterator end() const { return iterator(streams_.end()); }

    void push_back(ReadStream<ReadType> reader) {
        streams_.emplace_back(std::move(reader));
    }

    void reset() {
        for (auto &reader : streams_)
            reader.reset();
    }

    void close() {
        for (auto &reader : streams_)
            reader.close();
    }

    void clear() {
        streams_.clear();
    }
};
}
