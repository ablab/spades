/**
 * @file    scf_parser.hpp
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
 * ScfParser is the parser stream that reads data from .scf and .abi.
 * Be careful with this stream - it can read only one record from 
 * one file! It is experimental and is not suggected to be actively
 * used. 
 */

#ifndef COMMON_IO_SCFPARSER_HPP
#define COMMON_IO_SCFPARSER_HPP

#include <zlib.h>
#include <string>
#include <cassert>
#include <io_lib/Read.h>
#include "common/io/single_read.hpp"
#include "common/io/parser.hpp"
#include "common/sequence/quality.hpp"
#include "common/sequence/nucl.hpp"

namespace io {

class ScfParser : public Parser {
 public:
  /*
   * Default constructor.
   * 
   * @param filename The name of the file to be opened.
   * @param offset The offset of the read quality.
   */
  ScfParser(const std::string& filename,
         int offset = SingleRead::PHRED_OFFSET)
      :Parser(filename, offset), read_(NULL) {
    open();
  }

  /* 
   * Default destructor.
   */
  /* virtual */ ~ScfParser() {
    close();
  }

  /*
   * Read SingleRead from stream.
   *
   * @param read The SingleRead that will store read data.
   *
   * @return Reference to this stream.
   */
  /* virtual */ ScfParser& operator>>(SingleRead& read) {
    if (!is_open_ || eof_) {
      return *this;
    }
    // TODO(mariyafomkina): Rewrite 2 following instructions.
    read.SetName(filename_.c_str());     
    read.SetQuality("");
    read.SetSequence(read_->base);
    eof_ = true;
    return *this;
  }

  /*
   * Close the stream.
   */
  /* virtual */ void close() {
    if (is_open_) {
      is_open_ = false;
      eof_ = true;
    }
  }

 private:
  /*
   * @variable Data element that stores last SingleRead got from
   * stream.
   */ 
  Read* read_;

  /*
   * Open a stream.
   */
  /* virtual */ void open() {
    eof_ = false;
    is_open_ = true;
    read_ = read_reading(const_cast<char *>(filename_.c_str()), 0);
    if (read_ == NULLRead) {
      eof_ = true;
    }
  }

  /*
   * Hidden copy constructor.
   */
  ScfParser(const ScfParser& parser);
  /*
   * Hidden assign operator.
   */
  void operator=(const ScfParser& parser);
};

}

#endif /* COMMON_IO_SCFPARSER_HPP */
