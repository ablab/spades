////***************************************************************************
////* Copyright (c) 2011-2014 Saint-Petersburg Academic University
////* All Rights Reserved
////* See file LICENSE for details.
////****************************************************************************
//  todo remove!!!
///**
// * @file    cutting_reader_wrapper.hpp
// * @author  Mariya Fomkina
// * @version 1.0
// *
// * @section LICENSE
// *
// * This program is free software; you can redistribute it and/or
// * modify it under the terms of the GNU General Public License as
// * published by the Free Software Foundation; either version 2 of
// * the License, or (at your option) any later version.
// *
// * @section DESCRIPTION
// *
// * CuttingReaderWrapper is the class-wrapper that reads only set
// * number of reads from another reader.
// */
//
//#ifndef COMMON_IO_CUTTINGREADERWRAPPER_HPP_
//#define COMMON_IO_CUTTINGREADERWRAPPER_HPP_
//
//#include "io/ireader.hpp"
//
//namespace io {
//
//template<typename ReadType>
//class CuttingReaderWrapper : public IReader<ReadType> {
// public:
//  /*
//   * Default constructor.
//   *
//   * @param reader Reference to any other reader (child of IReader).
//   * @param cut Number of reads to be read (-1 by default, i.e. all).
//   */
//  explicit CuttingReaderWrapper(IReader<ReadType>& reader,
//                                size_t cut = -1)
//      : reader_(reader), cut_(cut), read_(0) {
//  }
//
//  /*
//   * Default destructor.
//   */
//  /* virtual */ ~CuttingReaderWrapper() {
//    close();
//  }
//
//  /*
//   * Check whether the stream is opened.
//   *
//   * @return true of the stream is opened and false otherwise.
//   */
//  /* virtual */ bool is_open() {
//    return reader_.is_open();
//  }
//
//  /*
//   * Check whether we've reached the end of stream.
//   *
//   * @return true if the end of stream is reached and false
//   * otherwise.
//   */
//  /* virtual */ bool eof() {
//    return (read_ == cut_) || (reader_.eof());
//  }
//
//  /*
//   * Read SingleRead or PairedRead from stream (according to ReadType).
//   *
//   * @param read The SingleRead or PairedRead that will store read
//   * data.
//   *
//   * @return Reference to this stream.
//   */
//  /* virtual */ CuttingReaderWrapper& operator>>(ReadType& read) {
//    if (read_ < cut_) {
//      reader_ >> read;
//      ++read_;
//    }
//    return (*this);
//  }
//
//  /*
//   * Close the stream.
//   */
//  /* virtual */ void close() {
//    reader_.close();
//  }
//
//  /*
//   * Close the stream and open it again.
//   */
//  /* virtual */ void reset() {
//    read_ = 0;
//    reader_.reset();
//  }
//
//  ReadStat get_stat() const {
//      return reader_.get_stat();
//  }
//
// private:
//  /*
//   * @variable Internal stream readers.
//   */
//  IReader<ReadType>& reader_;
//  /*
//   * @variable Number of reads that are allowed to read (if it is less
//   * than 0, all the reads in stream are allowed to be read).
//   */
//  size_t cut_;
//  /*
//   * @variable Number of reads that are read till the moment.
//   */
//  size_t read_;
//
//  /*
//   * Hidden copy constructor.
//   */
//  explicit CuttingReaderWrapper(const CuttingReaderWrapper<ReadType>&
//                                reader);
//  /*
//   * Hidden assign operator.
//   */
//  void operator=(const CuttingReaderWrapper<ReadType>& reader);
//};
//
//}
//
//#endif /* COMMON_IO_CUTTINGREADERWRAPPER_HPP_ */
