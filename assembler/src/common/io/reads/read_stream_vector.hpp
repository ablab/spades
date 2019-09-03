//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "read_stream.hpp"

#include <memory>
#include <vector>

namespace io {
//todo rename file

//todo check destroy_readers logic and usages
template<class ReadType>
class ReadStreamList {
public:
    typedef ReadType ReadT;
    typedef ReadStream<ReadType> ReaderT;

private:
    typedef std::vector<ReadStream<ReadType>> ReadersT;
    ReadersT readers_;

public:
    using iterator = typename ReadersT::iterator;
    using const_iterator = typename ReadersT::const_iterator;

    ReadStreamList() = default;
    ReadStreamList(const ReadStreamList&) = delete;
    ReadStreamList(ReadStreamList&&) = default;
    ReadStreamList &operator=(const ReadStreamList&) = delete;
    ReadStreamList &operator=(ReadStreamList&&) = default;
    
    explicit ReadStreamList(ReaderT reader_ptr) {
        readers_.emplace_back(std::move(reader_ptr));
    }

    ReadStream<ReadType> &operator[](size_t i) { return readers_[i]; }
    const ReadStream<ReadType> &operator[](size_t i) const { return readers_[i]; }

    ReaderT &back() {
        return readers_.back();
    }

    size_t size() const {
        return readers_.size();
    }

    bool eof() const {
        for (size_t i = 0; i < readers_.size(); ++i) {
            if (!readers_[i].eof()) {
                return false;
            }
        }
        return true;
    }

    iterator begin() { return readers_.begin(); }
    iterator end() { return iterator(readers_.end()); }
    const_iterator begin() const { return readers_.begin(); }
    const_iterator end() const { return iterator(readers_.end()); }

    void push_back(ReadStream<ReadType> reader) {
        readers_.emplace_back(std::move(reader));
    }

    void reset() {
        for (auto &reader : readers_)
            reader.reset();
    }

    void close() {
        for (auto &reader : readers_)
            reader.close();
    }

    void clear() {
        readers_.clear();
    }
};

}
