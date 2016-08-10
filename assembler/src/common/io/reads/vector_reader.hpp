//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once
#include "io/reads/ireadstream.hpp"
namespace io {

/**
 * Use vector<T> as input-stream with operator>>(T& t)
 */
template <typename T>
class VectorReadStream : public ReadStream<T> {
       std::vector<T> data_;
       size_t pos_;
       bool closed_;
public:
       VectorReadStream(const std::vector<T>& data) : data_(data), pos_(0), closed_(false) {

       }

       VectorReadStream(const T& item) : data_({item}), pos_(0), closed_(false) {

       }

       virtual bool eof() /*const */{
               return pos_ == data_.size();
       }

       VectorReadStream<T>& operator>>(T& t) {
           VERIFY(!eof());
           t = data_[pos_++];
           return *this;
       }

       void close() {
               closed_ = true;
       }

       virtual bool is_open() /*const */{
               return !closed_;
       }

       void reset() {
               pos_ = 0;
       }

       ReadStreamStat get_stat() const {
           //todo
           ReadStreamStat stat;
           stat.read_count_ = data_.size();

           return stat;
       }

};

}
