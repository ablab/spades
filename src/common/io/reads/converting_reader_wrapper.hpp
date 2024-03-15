//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "read_stream_vector.hpp"

namespace io {

/**
* SquashingWrapper is the class-wrapper that reads SingleReads
* from Reader<PairedRead> (first and second single reads in a pair
* one by one).
*/
template<class PairedReadType>
class SquashingWrapper {
    typedef typename PairedReadType::SingleReadT SingleReadT;
    typedef ReadStream<PairedReadType> PairedReader;
public:

    explicit SquashingWrapper(PairedReader reader)
            : reader_(std::move(reader)), pairedread_(), index_(0) {}

    /*
     * Check whether the stream is opened.
     *
     * @return true if the stream is opened and false otherwise.
     */
    bool is_open() {
        return reader_.is_open();
    }

    /*
     * Check whether we've reached the end of stream.
     *
     * @return true if the end of the stream is reached and false
     * otherwise.
     */
    bool eof() {
        return (index_ == 0) && (reader_.eof());
    }

    /*
     * Read SingleRead from stream (which is actually the part of
     * PairedRead from stream).
     *
     * @param singleread The SingleRead that will store read data.
     *
     * @return Reference to this stream.
     */
    SquashingWrapper &operator>>(SingleReadT &singleread) {
        if (index_ == 0)
            reader_ >> pairedread_;

        singleread = pairedread_[index_];
        index_ = 1 - index_;
        return *this;
    }

    /*
     * Close the stream.
     */
    void close() {
        reader_.close();
    }

    /*
     * Close the stream and open it again.
     */
    void reset() {
        index_ = 0;
        reader_.reset();
    }

    constexpr auto && unwrap() {
        return reader_;
    }
    template<class T>
    constexpr auto && recover() {
        return reader_.template recover<T>();
    }

private:
    /*
     * @variable Internal stream reader.
     */
    PairedReader reader_;
    /*
     * @variable Element that stores the last read PairedRead from
     * stream.
     */
    PairedReadType pairedread_;
    /*
     * @variable Index of current part of PairedRead.
     */
    size_t index_;

};

template<class PairedReadType>
ReadStream<typename PairedReadType::SingleReadT>
SquashingWrap(ReadStream<PairedReadType> reader_ptr) {
    return SquashingWrapper<PairedReadType>(std::move(reader_ptr));
}

template<class PairedReadType>
ReadStreamList<typename PairedReadType::SingleReadT> SquashingWrap(ReadStreamList<PairedReadType> readers) {
    ReadStreamList<typename PairedReadType::SingleReadT> answer;
    for (auto &reader : readers)
        answer.push_back(SquashingWrap<PairedReadType>(std::move(reader)));

    return answer;
}

}
