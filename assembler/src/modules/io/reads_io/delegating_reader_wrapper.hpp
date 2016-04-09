//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "ireader.hpp"

namespace io {

//todo rename file
template<typename ReadType>
class DelegatingWrapper: public ReadStream<ReadType> {
public:
    typedef std::shared_ptr<ReadStream<ReadType>> ReadStreamPtrT;

    explicit DelegatingWrapper(ReadStreamPtrT reader) : reader_(reader) {}


    /* virtual */ bool is_open() {
        return reader_->is_open();
    }

    /* virtual */ bool eof() {
        return reader_->eof();
    }

    /* virtual */ DelegatingWrapper& operator>>(ReadType& read) {
        (*reader_) >> read;
        return *this;
    }

    /* virtual */
    void close() {
        reader_->close();
    }

    /*
     * Close the stream and open it again.
     */
    /* virtual */
    void reset() {
        reader_->reset();
    }

    /* virtual */
    ReadStreamStat get_stat() const {
        return reader_->get_stat();
    }

protected:
    ReadStream<ReadType>& reader() {
        return *reader_;
    }

private:
    ReadStreamPtrT reader_;

};

}
