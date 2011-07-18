/**
 * @file    cutting_reader_wrapper.hpp
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
 * CuttingReaderWrapper is the class-wrapper that reads only set
 * number of reads from another reader.
 */

#ifndef COMMON_IO_CUTTINGREADERWRAPPER_HPP_
#define COMMON_IO_CUTTINGREADERWRAPPER_HPP_

#include "common/io/reader.hpp"

template<typename ReadType>
class CuttingReaderWrapper : public Reader<ReadType> {
 public:
  CuttingReaderWrapper(const Reader<ReadType>* reader,
                       size_t cut = -1)
      : reader_(reader), cut_(cut), read_(0) {
    is_open_ = (*reader_).is_open();
    eof_ = (*reader_).eof();
  }

  /*virtual*/~CuttingReaderWrapper() {
    close();
  }

  /*virtual*/CuttingReaderWrapper& operator>>(ReadType& read) {
    (*reader_) >> read;
    ++read_;
    if ((read_ == cut_) || ((*reader_).eof())) {
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
    read_ = 0;
    (*reader_).reset();
  }

 private:
  Reader<ReadType>* reader_;
  size_t read_;
  size_t cut_;
};

#endif /* COMMON_IO_CUTTINGREADERWRAPPER_HPP_ */
