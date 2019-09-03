//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "read_stream.hpp"

namespace io {

//todo rename file
template<typename ReadType>
class DelegatingWrapper {
public:
    typedef ReadStream<ReadType> ReadStreamT;

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

    constexpr auto && unwrap() {
        return reader_;
    }
    template<class T>
    constexpr auto && recover() {
        return reader_.template recover<T>();
    }

protected:
    ReadStream<ReadType>& reader() {
        return reader_;
    }

private:
    ReadStream<ReadType> reader_;
};

}
