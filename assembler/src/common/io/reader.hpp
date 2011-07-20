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
 * Reader is the base class that gets single reads or paired reads
 * from one or two input files respectively.
 * Reader<SingleRead> is the very base class that reads from one file
 * through Parser object.
 * Reader<PairedRead> is the class that reads data from two input
 * files and gets paired reads using this data and distance information.
 */

#ifndef COMMON_IO_READER_HPP_
#define COMMON_IO_READER_HPP_

#include "common/io/single_read.hpp"
#include "common/io/paired_read.hpp"
#include "common/io/parser.hpp"
#include "common/io/parser.cpp"  // TEMPORARY HACK!!!

/*
 * This class only represents Reader interface. All the functionality
 * is implemented in specializations. Thus, it's impossible to use
 * Reader<int> or Reader<std::string>. The only possible variants are
 * Reader<SingleRead> and Reader<PairedRead>.
 */
template<typename ReadType>
class Reader {
 public:
  Reader(const typename ReadType::FilenameType& filename,
         int offset = SingleRead::PHRED_OFFSET,
         size_t distance = 0) = 0;
  virtual ~Reader() = 0;
  virtual bool is_open() = 0;
  virtual bool eof() = 0;
  virtual Reader& operator>>(ReadType& read) = 0;
  virtual void close() = 0;
  virtual void reset() = 0;
};

template<>
class Reader<SingleRead> {
 public:
  Reader(const SingleRead::FilenameType& filename,
         int offset = SingleRead::PHRED_OFFSET,
         size_t distance = 0)
      : filename_(filename), offset_(offset) {
    parser_ = SelectParser(filename_, offset_);
  }

  virtual ~Reader() {
    close();
    delete parser_;
  }

  virtual bool is_open() {
    if (parser_ != NULL) {
      return parser_->is_open();
    } else {
      return false;
    }
  }

  virtual bool eof() {
    if (parser_ != NULL) {
      return parser_->eof();
    } else {
      return true;
    }
  }

  virtual Reader& operator>>(SingleRead& singleread) {
    if (parser_ != NULL) {
      (*parser_) >> singleread;
    }
    return *this;
  }

  virtual void close() {
    if (parser_ != NULL) {
      parser_->close();
    }
  }

  virtual void reset() {
    if (parser_ != NULL) {
      parser_->reset();
    }
  }

 private:
  SingleRead::FilenameType filename_;
  int offset_;
  Parser* parser_;

  explicit Reader(const Reader<SingleRead>& reader);
  void operator=(const Reader<SingleRead>& reader);
};

template<>
class Reader<PairedRead> {
 public:
  Reader(const PairedRead::FilenameType& filename,
         size_t distance,
         int offset = SingleRead::PHRED_OFFSET)
      : filename_(filename), distance_(distance), offset_(offset) {
    first_ = new Reader<SingleRead>(filename_.first, offset_);
    second_ = new Reader<SingleRead>(filename_.second, offset_);
  }

  virtual ~Reader() {
    close();
  }

  virtual bool is_open() {
    return first_->is_open() && second_->is_open();
  }

  virtual bool eof() {
    return first_->eof() || second_->eof();
  }

  virtual Reader& operator>>(PairedRead& pairedread) {
    SingleRead sr1, sr2;
    (*first_) >> sr1;
    (*second_) >> sr2;
    pairedread = PairedRead(sr1, sr2, distance_);  // is it correct?
    return *this;
  }

  virtual void close() {
    first_->close();
    second_->close();
  }

  virtual void reset() {
    first_->reset();
    second_->reset();
  }

 private:
  PairedRead::FilenameType filename_;
  size_t distance_;
  int offset_;
  Reader<SingleRead>* first_;
  Reader<SingleRead>* second_;

  explicit Reader(const Reader<PairedRead>& reader);
  void operator=(const Reader<PairedRead>& reader);
};

#endif /* COMMON_IO_READER_HPP_ */
