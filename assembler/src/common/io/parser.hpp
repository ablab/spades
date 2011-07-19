/**
 * @file    parser.hpp
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
 * Parser is the parent class for all streams that read data from
 * different file types (fastq, fasta, sam etc).
 */

#ifndef COMMON_IO_PARSER_HPP
#define COMMON_IO_PARSER_HPP

#include <string>
#include "common/io/single_read.hpp"

class Parser {
 public:
  Parser(const std::string& filename,
         int offset = SingleRead::PHRED_OFFSET)
      : filename_(filename), offset_(offset) {}

  virtual ~Parser() {}

  virtual bool is_open() const {
    return is_open_;
  }

  virtual bool eof() const {
    return eof_;
  }

  virtual Parser& operator>>(SingleRead& read) = 0;

  virtual void close() = 0;

  void reset() {
    close();
    open();
  }

 protected:
  std::string filename_;
  int offset_;
  bool is_open_;
  bool eof_;

 private:
  virtual void open() = 0;
};

Parser* SelectParser(const std::string filename_,
                     int offset = SingleRead::PHRED_OFFSET) {
  // get extension
  // select one of parsers according to extension
  // return parser*
}

#endif /* COMMON_IO_PARSER_HPP */
