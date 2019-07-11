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
class DelegatingWrapper {
public:
    typedef ReadStream<ReadType> ReadStreamPtrT;
    
    explicit DelegatingWrapper(ReadStream<ReadType> reader)
            : reader_{std::move(reader)} {}

    bool is_open() {
        return reader_.is_open();
    }

    bool eof() {
        return reader_.eof();
    }

    DelegatingWrapper& operator>>(ReadType& read) {
        reader_ >> read;
        return *this;
    }

    void close() {
        reader_.close();
    }

    void reset() {
        reader_.reset();
    }

protected:
    ReadStream<ReadType>& reader() {
        return reader_;
    }

private:
    ReadStream<ReadType> reader_;
};

}
