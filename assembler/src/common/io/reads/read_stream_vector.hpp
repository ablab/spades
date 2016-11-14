//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "ireader.hpp"
#include <vector>

namespace io {
//todo rename file

//todo check destroy_readers logic and usages
template<class ReadType>
class ReadStreamList {
public:
    typedef ReadType ReadT;
    typedef ReadStream<ReadType> ReaderT;
    typedef std::shared_ptr<ReaderT> ReaderPtrT;

private:
    std::vector<ReaderPtrT> readers_;

public:

    explicit ReadStreamList(const std::vector<ReaderPtrT> &readers) : readers_(readers) {
    }

    ReadStreamList() {
    }

    explicit ReadStreamList(ReaderT *reader_ptr) : readers_(1, ReaderPtrT(reader_ptr)) {
    }

    explicit ReadStreamList(ReaderPtrT reader_ptr) : readers_(1, reader_ptr) {
    }

    explicit ReadStreamList(size_t size) : readers_(size) {
    }

    //todo use boost iterator facade
    class iterator : public std::iterator<std::input_iterator_tag, ReaderT> {
        typedef typename std::vector<ReaderPtrT>::iterator vec_it;
        vec_it it_;
    public:

        iterator(vec_it it) : it_(it) {
        }

        void operator++() {
            ++it_;
        }

        bool operator==(const iterator &that) {
            return it_ == that.it_;
        }

        bool operator!=(const iterator &that) {
            return it_ != that.it_;
        }

        ReaderT &operator*() {
            return *(*it_);
        }
    };

    ReaderT &operator[](size_t i) {
        return *readers_.at(i);
    }

    ReaderPtrT &ptr_at(size_t i) {
        return readers_.at(i);
    }

    ReaderT &back() {
        return *readers_.back();
    }

    size_t size() const {
        return readers_.size();
    }

    bool eof() const {
        for (size_t i = 0; i < readers_.size(); ++i) {
            if (!readers_[i]->eof()) {
                return false;
            }
        }
        return true;
    }

    iterator begin() {
        return iterator(readers_.begin());
    }

    iterator end() {
        return iterator(readers_.end());
    }

    void push_back(ReaderT *reader_ptr) {
        readers_.push_back(ReaderPtrT(reader_ptr));
    }

    void push_back(ReaderPtrT reader_ptr) {
        readers_.push_back(reader_ptr);
    }

    void reset() {
        for (size_t i = 0; i < readers_.size(); ++i) {
            readers_[i]->reset();
        }
    }

    void close() {
        for (size_t i = 0; i < readers_.size(); ++i) {
            readers_[i]->close();
        }
    }

    void clear() {
        readers_.clear();
    }

    ReadStreamStat get_stat() const {
        ReadStreamStat stat;
        for (size_t i = 0; i < readers_.size(); ++i) {
            stat.merge(readers_[i]->get_stat());
        }
        return stat;
    }

};

}
