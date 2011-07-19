/**
 * @file    fastqgz_parser.hpp
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
 * FastqgzParser is the parser stream that reads data from .fastq.gz
 * files.
 */

#ifndef COMMON_IO_FASTQGZPARSER_HPP
#define COMMON_IO_FASTQGZPARSER_HPP

#include <string>
#include <cassert>
#include "libs/kseq/kseq.h"
#include <zlib.h>
#include "common/io/single_read.hpp"
#include "common/sequence/quality.hpp"
#include "common/sequence/nucl.hpp"

// STEP 1: declare the type of file handler and the read() function
KSEQ_INIT(gzFile, gzread)

class FastqgzParser : public Parser {
 public:
  FastqgzParser(const std::string& filename,
         int offset = SingleRead::PHRED_OFFSET)
      : Parser(filename, offset) {
    is_open_ = open(filename_);
  }

  /* virtual */ ~FastqgzParser() {
    close();
  }

  /* virtual */ FastqgzParser& operator>>(SingleRead& read) {
    // some strange reaction on the end of the stream
    // should be rewritten
    assert(is_open_);
    assert(!eof_);
    read.setName(seq_->name.s);
    if (seq_->qual.s) {
      read.setQuality(seq_->qual.s, offset_);
    }
    read.setSequence(seq_->seq.s);
    read_ahead(); // make actual read for the next result
    return *this;
  }

  /* virtual */ void close() {
    if (is_open_) {
      // STEP 5: destroy seq
      kseq_destroy(seq_);
      // STEP 6: close the file handler 
      gzclose(fp_); 
      is_open_ = false;
    }
  }

 private:
  gzFile fp_;
  kseq_t* seq_;

  /* virtual */ bool open() {
    // STEP 2: open the file handler
    fp_ = gzopen(filename_.c_str(), "r"); 
    if (!fp_) {
      return false;
    }
    is_open_ = true;
    // STEP 3: initialize seq
    seq_ = kseq_init(fp_); 
    eof_ = false;
    read_ahead();
    return true;
  }

  void read_ahead() {
    assert(is_open_);
    assert(!eof_);
    if (kseq_read(seq_) < 0) {
      eof_ = true;
    }
  }
};

#endif /* COMMON_IO_FASTQGZPARSER_HPP */
