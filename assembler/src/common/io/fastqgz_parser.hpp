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

#include <zlib.h>
#include <string>
#include <cassert>
#include "libs/kseq/kseq.h"
#include "common/io/single_read.hpp"
#include "common/io/parser.hpp"
#include "common/sequence/quality.hpp"
#include "common/sequence/nucl.hpp"

namespace fastqgz {
// STEP 1: declare the type of file handler and the read() function
KSEQ_INIT(gzFile, gzread)
}

class FastqgzParser : public Parser {
 public:
  FastqgzParser(const std::string& filename,
         int offset = SingleRead::PHRED_OFFSET)
      :Parser(filename, offset) {
    open();
  }

  /* virtual */ ~FastqgzParser() {
    close();
  }

  /* virtual */ FastqgzParser& operator>>(SingleRead& read) {
    // some strange reaction on the end of the stream
    // should be rewritten
    assert(is_open_);
    assert(!eof_);
    read.SetName(seq_->name.s);
    if (seq_->qual.s) {
      read.SetQuality(seq_->qual.s, offset_);
    }
    read.SetSequence(seq_->seq.s);
    ReadAhead();
    return *this;
  }

  /* virtual */ void close() {
    if (is_open_) {
      // STEP 5: destroy seq
      fastqgz::kseq_destroy(seq_);
      // STEP 6: close the file handler
      gzclose(fp_);
      is_open_ = false;
    }
  }

 private:
  gzFile fp_;
  fastqgz::kseq_t* seq_;

  /* virtual */ void open() {
    // STEP 2: open the file handler
    fp_ = gzopen(filename_.c_str(), "r");
    if (!fp_) {
      is_open_ = false;
      return;
    }
    // STEP 3: initialize seq
    seq_ = fastqgz::kseq_init(fp_);
    eof_ = false;
    is_open_ = true;
    ReadAhead();
  }

  void ReadAhead() {
    assert(is_open_);
    assert(!eof_);
    if (fastqgz::kseq_read(seq_) < 0) {
      eof_ = true;
    }
  }
};

#endif /* COMMON_IO_FASTQGZPARSER_HPP */
