/**
 * @file    sff_parser.hpp
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
 * SamBamParser is the parser stream that reads data from .sam and
 * .bam files.
 */

#ifndef COMMON_IO_SFFPARSER_HPP
#define COMMON_IO_SFFPARSER_HPP

#include <zlib.h>
#include <string>
#include <cassert>
#include <io_lib/Read.h>
#include "common/io/single_read.hpp"
#include "common/io/parser.hpp"
#include "common/sequence/quality.hpp"
#include "common/sequence/nucl.hpp"

namespace io {

class SffParser : public Parser {
 public:
  /*
   * Default constructor.
   * 
   * @param filename The name of the file to be opened.
   * @param offset The offset of the read quality.
   */
  SffParser(const std::string& filename,
         int offset = SingleRead::PHRED_OFFSET)
      :Parser(filename, offset), read_(NULL) {
    open();
  }

  /* 
   * Default destructor.
   */
  /* virtual */ ~SffParser() {
    close();
  }

  /*
   * Read SingleRead from stream.
   *
   * @param read The SingleRead that will store read data.
   *
   * @return Reference to this stream.
   */
  /* virtual */ SffParser& operator>>(SingleRead& read) {
    if (!is_open_ || eof_) {
      return *this;
    }
    // read.SetName(seq_.getReadName());
    // if (seq_.getQuality()) {
    //   read.SetQuality(seq_.getQuality(), offset_);
    // }
    // read.SetSequence(seq_.getSequence());
    read.SetName("Name");
    read.SetQuality("");
    read.SetSequence(read_->base);
    ReadAhead();
    return *this;
  }

  /*
   * Close the stream.
   */
  /* virtual */ void close() {
    if (is_open_) {
      //fp_.Close();
      is_open_ = false;
      eof_ = true;
    }
  }

 private:
  /* 
   * @variable File that is associated with gzipped data file.
   */
  //SamFile fp_;
  /*
   * @variable File header.
   */
  //SamFileHeader header_;
  /*
   * @variable Data element that stores last SingleRead got from
   * stream.
   */ 
  //SamRecord seq_;
  Read* read_;

  /*
   * Open a stream.
   */
  /* virtual */ void open() {
    //if (fp_.OpenForRead(filename_.c_str())) {
    //  fp_.ReadHeader(header_);
      eof_ = false;
      is_open_ = true;
      ReadAhead();
      //}
  }

  /* 
   * Read next SingleRead from file.
   */
  void ReadAhead() {
    assert(is_open_);
    assert(!eof_);
    //if (fp_.ReadRecord(header_, seq_) == 0) {
    //  eof_ = true;
    //}
    read_ = read_reading(const_cast<char *>(filename_.c_str()), 0);
    if (read_ == NULLRead) {
      eof_ = true;
    }
  }

  /*
   * Hidden copy constructor.
   */
  SffParser(const SffParser& parser);
  /*
   * Hidden assign operator.
   */
  void operator=(const SffParser& parser);
};

}

#endif /* COMMON_IO_SFFPARSER_HPP */
