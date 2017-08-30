//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <vector>
#include <string>
#include <cstdlib>
#include <cstdint>

namespace utils {

template<class V>
class ValueArray {
    static const size_t InvalidIdx = SIZE_MAX;
public:
    typedef size_t IdxType;
    typedef V ValueType;

protected:
    typedef std::vector<V> StorageT;
    StorageT data_;

    void resize(size_t size) {
        data_.resize(size);
    }

public:
    typedef typename StorageT::iterator value_iterator;
    typedef typename StorageT::const_iterator const_value_iterator;

    ValueArray() {
    }

    ~ValueArray() {
    }

    void clear() {
        data_.clear();
        StorageT().swap(data_);
    }

    const V &operator[](size_t idx) const {
        return data_[idx];
    }

    V &operator[](size_t idx) {
        return data_[idx];
    }

public:
    size_t size() const {
        return data_.size();
    }

    value_iterator value_begin() {
        return data_.begin();
    }
    const_value_iterator value_begin() const {
        return data_.begin();
    }
    const_value_iterator value_cbegin() const {
        return data_.cbegin();
    }
    value_iterator value_end() {
        return data_.end();
    }
    const_value_iterator value_end() const {
        return data_.end();
    }
    const_value_iterator value_cend() const {
        return data_.cend();
    }

    template<class Writer>
    void BinWrite(Writer &writer) const {
        size_t sz = data_.size();
        writer.write((char*) &sz, sizeof(sz));
        writer.write((char*) &data_[0], sz * sizeof(data_[0]));
    }

    template<class Reader>
    void BinRead(Reader &reader, const std::string &) {
        clear();
        size_t sz = 0;
        reader.read((char*) &sz, sizeof(sz));
        data_.resize(sz);
        reader.read((char*) &data_[0], sz * sizeof(data_[0]));
    }
};

}
