//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

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
 * ConvertingReaderWrapper is the class-wrapper that reads SingleReads 
 * from IReader<PairedRead> (first and second single reads in a pair
 * one by one).
 */

#ifndef COMMON_IO_CONVERTINGREADERWRAPPER_HPP_
#define COMMON_IO_CONVERTINGREADERWRAPPER_HPP_

#include "io/single_read.hpp"
#include "io/paired_read.hpp"
#include "io/ireader.hpp"

namespace io {

class ConvertingReaderWrapper : public IReader<SingleRead> {
 public:
  /*
   * Default constructor.
   *
   * @param reader Reference to any other reader (child of IReader).
   */
  explicit ConvertingReaderWrapper(IReader<PairedRead>& reader)
      : reader_(reader), pairedread_(), index_(0) {
  }

  /* 
   * Default destructor.
   */
  /* virtual */~ConvertingReaderWrapper() {
    close();
  }

  /* 
   * Check whether the stream is opened.
   *
   * @return true if the stream is opened and false otherwise.
   */
  /* virtual */ bool is_open() {
    return reader_.is_open();
  }

  /* 
   * Check whether we've reached the end of stream.
   *
   * @return true if the end of the stream is reached and false
   * otherwise.
   */
  /* virtual */ bool eof() {
    return (index_ == 0) && (reader_.eof());
  }

  /*
   * Read SingleRead from stream (which is actually the part of
   * PairedRead from stream).
   *
   * @param singleread The SingleRead that will store read data.
   *
   * @return Reference to this stream.
   */
  /* virtual */ ConvertingReaderWrapper& operator>>(
      SingleRead& singleread) {
    if (index_ == 0) {
      reader_ >> pairedread_;
    }
    singleread = pairedread_[index_];
    index_ = 1 - index_;
    return (*this);
  }

  /*
   * Close the stream.
   */
  /* virtual */ void close() {
    reader_.close();
  }

  /* 
   * Close the stream and open it again.
   */
  /* virtual */ void reset() {
    index_ = 0;
    reader_.reset();
  }

  ReadStat get_stat() const {
        return reader_.get_stat();
  }

 private:
  /*
   * @variable Internal stream reader.
   */
  IReader<PairedRead>& reader_;
  /*
   * @variable Element that stores the last read PairedRead from
   * stream.
   */
  PairedRead pairedread_;
  /*
   * @variable Index of current part of PairedRead.
   */
  size_t index_;

  /*
   * Hidden copy constructor.
   */
  explicit ConvertingReaderWrapper(const ConvertingReaderWrapper&
                                   reader);
  /*
   * Hidden assign operator.
   */
  void operator=(const ConvertingReaderWrapper& reader);
};

}

#endif /* COMMON_IO_CONVERTINGREADERWRAPPER_HPP_ */
