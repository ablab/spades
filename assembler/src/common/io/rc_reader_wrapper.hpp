/**
 * @file    rc_reader_wrapper.hpp
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
 * RCReaderWrapper is the class-wrapper that gets reads and reverse
 * complimentary reads from given reader (one by one).
 */

#ifndef COMMON_IO_RCREADERWRAPPER_HPP_
#define COMMON_IO_RCREADERWRAPPER_HPP_

#include "common/io/reader.hpp"

template<typename ReadType>
class RCReaderWrapper : public Reader<ReadType> {
 public:
  explicit RCReaderWrapper(const Reader<ReadType>* reader)
      : reader_(reader), was_rc_(true) {
    is_open_ = (*reader_).is_open();
    eof_ = (*reader_).eof();
  }

  /*virtual*/~RCReaderWrapper() {
    close();
  }

  /*virtual*/RCReaderWrapper& operator>>(ReadType& read) {
    if (was_rc_) {
      (*reader_) >> read;
      rc_read_ = read;
    } else {
      read = !rc_read_;
    }
    was_rc_ = !was_rc_;
    if ((was_rc_) && ((*reader_).eof())) {
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
    was_rc_ = false;
    (*reader_).reset();
  }

 private:
  Reader<ReadType>* reader_;
  ReadType rc_read_;
  bool was_rc_;
};

#endif /* COMMON_IO_RCREADERWRAPPER_HPP_ */
