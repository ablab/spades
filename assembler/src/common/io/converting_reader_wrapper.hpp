/**
 * @file    converting_reader_wrapper.hpp
 * @author  Mariya Fomkina
 * @version 1.0
 *
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 *
 * @section DESCRIPTION
 *
 * ConvertingReaderWrapper is the class-wrapper that reads single
 * reads from paired reader (first and second single reads in a pair
 * one by one.
 */

#ifndef COMMON_IO_CONVERTINGREADERWRAPPER_HPP_
#define COMMON_IO_CONVERTINGREADERWRAPPER_HPP_

#include "common/io/reader.hpp"

class ConvertingReaderWrapper : public Reader<SingleRead> {
 public:
  explicit ConvertingReaderWrapper(const Reader<PairedRead>* reader)
      : reader_(reader), index_(0) {
    is_open_ = (*reader_).is_open();
    eof_ = (*reader_).eof();
  }

  /*virtual*/~ConvertingReaderWrapper() {
    close();
  }

  /*virtual*/ConvertingReaderWrapper& operator>>(SingleRead& singleread) {
    if (index_ == 0) {
      reader_ >> pairedread_;
    }
    singleread = pairedread_[index_];
    index_ = 1 - index_;
    if ((*reader_).eof()) {
      eof_ = true;
    }
    return (*this);
  }

  /*virtual*/void close() {
    if (is_open_) {
      (*reader_).close();
      is_open_ = false;
    }
  }

  /*virtual*/void reset() {
    index_ = 0;
    (*reader_).reset();
  }

 private:
  Reader<PairedRead>* reader_;
  PairedRead pairedread_;
  size_t index_;
};

#endif /* COMMON_IO_CONVERTINGREADERWRAPPER_HPP_ */
