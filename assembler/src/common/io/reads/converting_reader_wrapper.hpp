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
class SquashingWrapper : public ReadStream<typename PairedReadType::SingleReadT> {
    typedef typename PairedReadType::SingleReadT SingleReadT;
    typedef std::shared_ptr<ReadStream<PairedReadType>> PairedReaderPtrT;
public:

    explicit SquashingWrapper(PairedReaderPtrT reader)
            : reader_(reader), pairedread_(), index_(0) {
    }

    /*
     * Check whether the stream is opened.
     *
     * @return true if the stream is opened and false otherwise.
     */
    /* virtual */ bool is_open() {
        return reader_->is_open();
    }

    /*
     * Check whether we've reached the end of stream.
     *
     * @return true if the end of the stream is reached and false
     * otherwise.
     */
    /* virtual */ bool eof() {
        return (index_ == 0) && (reader_->eof());
    }

    /*
     * Read SingleRead from stream (which is actually the part of
     * PairedRead from stream).
     *
     * @param singleread The SingleRead that will store read data.
     *
     * @return Reference to this stream.
     */
    /* virtual */ SquashingWrapper &operator>>(
            SingleReadT &singleread) {
        if (index_ == 0) {
            (*reader_) >> pairedread_;
        }
        singleread = pairedread_[index_];
        index_ = 1 - index_;
        return (*this);
    }

    /*
     * Close the stream.
     */
    /* virtual */ void close() {
        reader_->close();
    }

    /*
     * Close the stream and open it again.
     */
    /* virtual */ void reset() {
        index_ = 0;
        reader_->reset();
    }

private:
    /*
     * @variable Internal stream reader.
     */
    PairedReaderPtrT reader_;
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
std::shared_ptr<ReadStream<typename PairedReadType::SingleReadT>> SquashingWrap(
        std::shared_ptr<ReadStream<PairedReadType>> reader_ptr) {
    return std::make_shared<SquashingWrapper<PairedReadType>>(reader_ptr);
}

template<class PairedReadType>
ReadStreamList<typename PairedReadType::SingleReadT> SquashingWrap(const ReadStreamList<PairedReadType> &readers) {
    ReadStreamList<typename PairedReadType::SingleReadT> answer;
    for (size_t i = 0; i < readers.size(); ++i) {
        answer.push_back(SquashingWrap<PairedReadType>(readers.ptr_at(i)));
    }
    return answer;
}

//template<class ReaderPtrType>
//std::shared_ptr<Reader<typename ReaderPtrType::element_type::ReadT::SingleReadT>> SquashingWrap(ReaderPtrType reader_ptr) {
//    return std::make_shared<SquashingWrapper<typename ReaderPtrType::element_type::ReadT>>(reader_ptr);
//}
}
