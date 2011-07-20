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
  /*
   * Default constructor.
   * 
   * @param filename The name of the file to be opened.
   * @param offset The offset of the read quality.
   */
  Parser(const std::string& filename,
         int offset = SingleRead::PHRED_OFFSET)
      : filename_(filename), offset_(offset) {}

  /* 
   * Default destructor.
   */
  virtual ~Parser() {}

  /* 
   * Check whether the stream is opened.
   */
  virtual bool is_open() const {
    return is_open_;
  }

  /* 
   * Check whether we've reached the end of stream.
   */
  virtual bool eof() const {
    return eof_;
  }

  /*
   * Read single read from stream.
   *
   * @param read The single read that will store read data.
   *
   * @return Reference to this stream.
   */
  virtual Parser& operator>>(SingleRead& read) = 0;

  /*
   * Close the stream.
   */
  virtual void close() = 0;

  /* 
   * Close stream and open it again.
   */
  void reset() {
    close();
    open();
  }

 protected:
  /* 
   * @variable The name the file which stream reads from.
   */
  std::string filename_;
  /*
   * @variable Quality offset.
   */
  int offset_;
  /*
   * @variable Flag that shows whether the stream is opened.
   */
  bool is_open_;
  /*
   * @variable Flag that shows whether the end of the stream is
   * reached.
   */
  bool eof_;

 private:
  /* 
   * Open a stream.
   */
  virtual void open() = 0;
};

/*
 * Get extension from filename.
 *
 * @param filename The name of the file to read from.
 *
 * @return File extension (e.g. "fastq", "fastq.gz").
 */
std::string GetExtension(const std::string& filename);

/*
 * Select parser type according to file extension.
 *
 * @param filename The name of the file to be opened.
 * @param offset The offset of the read quality.

 * @return Pointer to the new parser object with these filename and
 * offset.
 */
Parser* SelectParser(const std::string& filename,
                     int offset = SingleRead::PHRED_OFFSET);

#endif /* COMMON_IO_PARSER_HPP */
