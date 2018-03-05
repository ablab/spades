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

    explicit DelegatingWrapper(
            const DelegatingWrapper& reader) = delete;

    void operator=(const DelegatingWrapper& reader) = delete;

    bool is_open() override {
        return reader_->is_open();
    }

    bool eof() override {
        return reader_->eof();
    }

    DelegatingWrapper& operator>>(ReadType& read) override {
        (*reader_) >> read;
        return *this;
    }

    void close() override {
        reader_->close();
    }

    /*
     * Close the stream and open it again.
     */
    void reset() override {
        reader_->reset();
    }

protected:
    ReadStream<ReadType>& reader() {
        return *reader_;
    }

private:
    ReadStreamPtrT reader_;
};

}
