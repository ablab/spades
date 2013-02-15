//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/**
 * @file    reader.hpp
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
 * Reader is the base class that gets SingleReads or PairedReads
 * from one or two input files respectively.
 * Reader<SingleRead> is the very base class that reads from one file
 * through Parser object.
 * Reader<PairedRead> is the class that reads data from two input
 * files and gets PairedReads using this data and distance information.
 */

#ifndef COMMON_IO_READER_HPP_
#define COMMON_IO_READER_HPP_

#include "io/ireader.hpp"
#include "io/single_read.hpp"
#include "io/paired_read.hpp"
#include "io/parser.hpp"
#include "simple_tools.hpp"

namespace io {

class Reader : public IReader<SingleRead> {
 public:
  /*
   * Default constructor.
   * 
   * @param filename The name of the file to be opened.
   * @param distance Doesn't have any sense here, but necessary for
   * wrappers.
   * @param offset The offset of the read quality.
   */
  explicit Reader(const SingleRead::FilenameType& filename,
                  OffsetType offset_type = PhredOffset)
      : filename_(filename), offset_type_(offset_type), parser_(NULL) {
	  CheckFileExistenceFATAL(filename_);
	  parser_ = SelectParser(filename_, offset_type_);
  }

  /* 
   * Default destructor.
   */
  /* virtual */ ~Reader() {
    close();
    delete parser_;
  }

  /* 
   * Check whether the stream is opened.
   *
   * @return true of the stream is opened and false otherwise.
   */
  /* virtual */ bool is_open() {
    if (parser_ != NULL) {
      return parser_->is_open();
    } else {
      return false;
    }
  }

  /* 
   * Check whether we've reached the end of stream.
   *
   * @return true if the end of stream is reached and false
   * otherwise.
   */
  /* virtual */ bool eof() {
    if (parser_ != NULL) {
      return parser_->eof();
    } else {
      return true;
    }
  }

  /*
   * Read SingleRead from stream.
   *
   * @param singleread The SingleRead that will store read data.
   *
   * @return Reference to this stream.
   */
  /* virtual */ Reader& operator>>(SingleRead& singleread) {
    if (parser_ != NULL) {
      (*parser_) >> singleread;
    }
    return *this;
  }

  /*
   * Close the stream.
   */
  /* virtual */ void close() {
    if (parser_ != NULL) {
      parser_->close();
    }
  }

  /* 
   * Close the stream and open it again.
   */
  /* virtual */ void reset() {
    if (parser_ != NULL) {
      parser_->reset();
    }
  }

  ReadStat get_stat() const {
        return ReadStat();
  }

 private:
  /* 
   * @variable The name of the file which stream reads from.
   */
  SingleRead::FilenameType filename_;
  /*
   * @variable Quality offset type.
   */
  OffsetType offset_type_;
  /*
   * @variable Internal stream that reads from file.
   */ 
  Parser* parser_;

  /*
   * Hidden copy constructor.
   */
  explicit Reader(const Reader& reader);
  /*
   * Hidden assign operator.
   */
  void operator=(const Reader& reader);


};


class SeparateReader : public IReader<PairedRead> {
 public:
  /*
   * Default constructor.
   *
   * @param filename The pair that containes the names of two files to
   * be opened.
   * @param distance Distance between parts of PairedReads.
   * @param offset The offset of the read quality.
   */
  explicit SeparateReader(const PairedRead::FilenamesType& filenames,
         size_t insert_size, bool change_order = false,
         bool revert_second = true,
         OffsetType offset_type = PhredOffset)
      : filenames_(filenames), insert_size_(insert_size),
        change_order_(change_order),
        revert_second_(revert_second),
        offset_type_(offset_type),
        first_(new Reader(filenames_.first, offset_type_)),
        second_(new Reader(filenames_.second, offset_type_)) {}

  /*
   * Default destructor.
   */
  /* virtual */ ~SeparateReader() {
    close();
    delete first_;
    delete second_;
  }

  /*
   * Check whether the stream is opened.
   *
   * @return true of the stream is opened and false otherwise.
   */
  /* virtual */ bool is_open() {
    return first_->is_open() && second_->is_open();
  }

  /*
   * Check whether we've reached the end of stream.
   *
   * @return true if the end of stream is reached and false
   * otherwise.
   */
  /* virtual */ bool eof() {
    return first_->eof() || second_->eof();
  }

  /*
   * Read PairedRead from stream.
   *
   * @param pairedread The PairedRead that will store read data.
   *
   * @return Reference to this stream.
   */
  /* virtual */ SeparateReader& operator>>(PairedRead& pairedread) {
    SingleRead sr1, sr2;
    (*first_) >> sr1;
    (*second_) >> sr2;
    if (revert_second_) sr2 = !sr2;

    pairedread = change_order_ ? PairedRead(sr2, sr1, insert_size_) : PairedRead(sr1, sr2, insert_size_);
    return *this;
  }

  /*
   * Close the stream.
   */
  /* virtual */ void close() {
    first_->close();
    second_->close();
  }

  /*
   * Close the stream and open it again.
   */
  /* virtual */ void reset() {
    first_->reset();
    second_->reset();
  }

  ReadStat get_stat() const {
    return ReadStat();
  }

 private:
  /*
   * @variable The names of the files which stream reads from.
   */
  PairedRead::FilenamesType filenames_;

  size_t insert_size_;

  bool change_order_;

  bool revert_second_;

  /*
   * @variable Quality offset type.
   */
  OffsetType offset_type_;

  /*
   * @variable The first stream (reads from first file).
   */
  Reader* first_;
  /*
   * @variable The second stream (reads from second file).
   */
  Reader* second_;

  /*
   * Hidden copy constructor.
   */
  explicit SeparateReader(const SeparateReader& reader);
  /*
   * Hidden assign operator.
   */
  void operator=(const SeparateReader& reader);


};


class MixedReader : public IReader<PairedRead> {
 public:
  /*
   * Default constructor.
   *
   * @param filename Single file
   * @param distance Distance between parts of PairedReads.
   * @param offset The offset of the read quality.
   */
  explicit MixedReader(const std::string& filename, size_t insert_size, bool change_order = false,
		  bool revert_second = true, OffsetType offset_type = PhredOffset)
      : filename_(filename), insert_size_(insert_size),
        change_order_(change_order),
        revert_second_(revert_second), offset_type_(offset_type),
        single_(new Reader(filename_, offset_type_)) {}

  /*
   * Default destructor.
   */
  /* virtual */ ~MixedReader() {
    close();
    delete single_;
  }

  /*
   * Check whether the stream is opened.
   *
   * @return true of the stream is opened and false otherwise.
   */
  /* virtual */ bool is_open() {
    return single_->is_open();
  }

  /*
   * Check whether we've reached the end of stream.
   *
   * @return true if the end of stream is reached and false
   * otherwise.
   */
  /* virtual */ bool eof() {
    return single_->eof();
  }

  /*
   * Read PairedRead from stream.
   *
   * @param pairedread The PairedRead that will store read data.
   *
   * @return Reference to this stream.
   */
  /* virtual */ MixedReader& operator>>(PairedRead& pairedread) {
    SingleRead sr1, sr2;
    (*single_) >> sr1;
    (*single_) >> sr2;

    if (revert_second_) sr2 = !sr2;

    pairedread = change_order_ ? PairedRead(sr2, sr1, insert_size_) : PairedRead(sr1, sr2, insert_size_);

    return *this;
  }

  /*
   * Close the stream.
   */
  /* virtual */ void close() {
    single_->close();
  }

  /*
   * Close the stream and open it again.
   */
  /* virtual */ void reset() {
    single_->reset();
  }

  ReadStat get_stat() const {
        return ReadStat();
  }

 private:
  /*
   * @variable The names of the file which stream reads from.
   */
  std::string filename_;

  size_t insert_size_;

  bool change_order_;

  bool revert_second_;

  /*
   * @variable Quality offset type.
   */
  OffsetType offset_type_;

  /*
   * @variable The single read stream.
   */
  Reader* single_;

  /*
   * Hidden copy constructor.
   */
  explicit MixedReader(const MixedReader& reader);
  /*
   * Hidden assign operator.
   */
  void operator=(const MixedReader& reader);

};


}

#endif /* COMMON_IO_READER_HPP_ */
