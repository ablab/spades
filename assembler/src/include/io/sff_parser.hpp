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
 * SffParser is a parser that reads data from .sff files.
 * Note: Quality offset type is set to PhredOffset.
 * Note: Current version of parser doesn't trim ends!
 */

#ifndef COMMON_IO_SFFPARSER_HPP
#define COMMON_IO_SFFPARSER_HPP

#include <zlib.h>
#include <string>
#include <cassert>
#include <io_lib/sff.h>
#include "io/single_read.hpp"
#include "io/parser.hpp"
#include "sequence/quality.hpp"
#include "sequence/nucl.hpp"

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
            OffsetType offset_type = PhredOffset)
      :Parser(filename, offset_type), read_() {
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
    read.SetName(read_.name().c_str());
    read.SetQuality(read_.GetPhredQualityString().c_str());
    read.SetSequence(read_.GetSequenceString().c_str());
    if (cnt_ == num_of_reads_) {
      eof_ = true;
    } else {
      ReadAhead();
    }    
    return *this;
  }

  /*
   * Close the stream.
   */
  /* virtual */ void close() {
    if (is_open_) {
      mfclose(sff_fp_);
      free_sff_common_header(h_);
      is_open_ = false;
      eof_ = true;
    }
  }

 private:
  /*
   * Header of the whole sff file.
   */
  sff_common_header* h_;
  /*
   * Header of one read in sff file.
   */
  sff_read_header* rh_;
  /*
   * Read data of one read in sff file.
   */
  sff_read_data* rd_;
  /*
   * File that contains sff data.
   */
  mFILE *sff_fp_;
  /* 
   * Number of reads that are stored in this sff file.
   */
  int num_of_reads_;
  /*
   * Number of reads that are already read from stream.
   */
  int cnt_;
  /*
   * Read that was last read from file (and next to be outputted by
   * operator>>.
   */
  SingleRead read_;

  /*
   * Open a stream.
   */
  /* virtual */ void open() {
    sff_fp_ = mfopen(filename_.c_str(), "r");
    if (sff_fp_ == NULL) {
      eof_ = true;
      is_open_ = false;
    } else {
      h_ = read_sff_common_header(sff_fp_);
      num_of_reads_ = (int) h_->nreads;
      cnt_ = 0;
      eof_ = false;
      is_open_ = true;
      ReadAhead();
    }
  }

  /* 
   * Read next SingleRead from file.
   */
  void ReadAhead() {
    assert(is_open_);
    assert(!eof_);
    rh_ = read_sff_read_header(sff_fp_);
    rd_ = read_sff_read_data(sff_fp_, h_->flow_len, rh_->nbases);
    ++cnt_;

    char* name_;
    char* bases_;
    unsigned char* quality_;
    unsigned char quality_char_;
    size_t nlen = rh_->name_len;
    size_t slen = rh_->nbases;
    name_ = (char *)malloc(nlen + 1);
    strncpy(name_, rh_->name, nlen);
    name_[nlen] = '\0';
    read_.SetName(name_);
    bases_ = (char *)malloc(slen + 1);
    strncpy(bases_, rd_->bases, slen);
    bases_[slen] = '\0';
    read_.SetSequence(bases_);
    quality_ = (unsigned char *)malloc(slen + 1);
    for (size_t i = 0; i < slen; i++) {
      quality_char_ = (rd_->quality[i] <= 93 ? rd_->quality[i] : 93) + 33;
      quality_[i] = quality_char_;
    }
    quality_[slen] = '\0';
    read_.SetQuality((char *)quality_, PhredOffset);
    free(name_);
    free(bases_);
    free(quality_);
    free_sff_read_header(rh_);
    free_sff_read_data(rd_);
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
