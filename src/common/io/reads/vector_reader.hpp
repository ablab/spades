//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
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
class VectorReadStream  {
       std::vector<T> data_;
       size_t pos_;
       bool closed_;
public:
       VectorReadStream(const std::vector<T>& data)
                     : data_(data), pos_(0), closed_(false) {}
       
       VectorReadStream(const T& item)
                     : data_({item}), pos_(0), closed_(false) {}

       bool eof() /*const */{
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

       bool is_open() /*const */{
              return !closed_;
       }

       void reset() {
               pos_ = 0;
       }

};

}
