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
#include "io/single_read.hpp"
#include "io/fasta_fastq_gz_parser.hpp"

namespace io {

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
      : filename_(filename), offset_(offset),
        is_open_(false), eof_(true) {}

  /* 
   * Default destructor.
   */
  virtual ~Parser() {}

  /* 
   * Check whether the stream is opened.
   *
   * @return true of the stream is opened and false otherwise.
   */
  virtual bool is_open() const {
    return is_open_;
  }

  /* 
   * Check whether we've reached the end of stream.
   *
   * @return true if the end of stream is reached and false
   * otherwise.
   */
  virtual bool eof() const {
    return eof_;
  }

  /*
   * Read SingleRead from stream.
   *
   * @param read The SingleRead that will store read data.
   *
   * @return Reference to this stream.
   */
  virtual Parser& operator>>(SingleRead& read) = 0;

  /*
   * Close the stream.
   */
  virtual void close() = 0;

  /* 
   * Close the stream and open it again.
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

class FastaFastqGzParser : public Parser {
 public:
  /*
   * Default constructor.
   *
   * @param filename The name of the file to be opened.
   * @param offset The offset of the read quality.
   */
  FastaFastqGzParser(const std::string& filename,
         int offset = SingleRead::PHRED_OFFSET)
      :Parser(filename, offset), fp_(), seq_(NULL) {
    open();
  }

  /*
   * Default destructor.
   */
  /* virtual */ ~FastaFastqGzParser() {
    close();
  }

  /*
   * Read SingleRead from stream.
   *
   * @param read The SingleRead that will store read data.
   *
   * @return Reference to this stream.
   */
  /* virtual */ FastaFastqGzParser& operator>>(SingleRead& read) {
    if (!is_open_ || eof_) {
      return *this;
    }
    read.SetName(seq_->name.s);
    if (seq_->qual.s) {
      read.SetQuality(seq_->qual.s, offset_);
    }
    read.SetSequence(seq_->seq.s);
    ReadAhead();
    return *this;
  }

  /*
   * Close the stream.
   */
  /* virtual */ void close() {
    if (is_open_) {
      // STEP 5: destroy seq
      fastafastqgz::kseq_destroy(seq_);
      // STEP 6: close the file handler
      gzclose(fp_);
      is_open_ = false;
      eof_ = true;
    }
  }

 private:
  /*
   * @variable File that is associated with gzipped data file.
   */
  gzFile fp_;
  /*
   * @variable Data element that stores last SingleRead got from
   * stream.
   */
  fastafastqgz::kseq_t* seq_;

  /*
   * Open a stream.
   */
  /* virtual */ void open() {
    // STEP 2: open the file handler
    fp_ = gzopen(filename_.c_str(), "r");
    if (!fp_) {
      is_open_ = false;
      return;
    }
    // STEP 3: initialize seq
    seq_ = fastafastqgz::kseq_init(fp_);
    eof_ = false;
    is_open_ = true;
    ReadAhead();
  }

  /*
   * Read next SingleRead from file.
   */
  void ReadAhead() {
    assert(is_open_);
    assert(!eof_);
    if (fastafastqgz::kseq_read(seq_) < 0) {
      eof_ = true;
    }
  }

  /*
   * Hidden copy constructor.
   */
  FastaFastqGzParser(const FastaFastqGzParser& parser);
  /*
   * Hidden assign operator.
   */
  void operator=(const FastaFastqGzParser& parser);
};

/*
 * Get extension from filename.
 *
 * @param filename The name of the file to read from.
 *
 * @return File extension (e.g. "fastq", "fastq.gz").
 */
inline std::string GetExtension(const std::string& filename) {
	  std::string name = filename;
	  size_t pos = name.find_last_of(".");
	  std::string ext = "";
	  if (pos != std::string::npos) {
	    ext = name.substr(name.find_last_of(".") + 1);
	    if (ext == "gz") {
	      ext = name.substr(name.find_last_of
	                        (".", name.find_last_of(".") - 1) + 1);
	    }
	  }
	  return ext;
}

/*
 * Select parser type according to file extension.
 *
 * @param filename The name of the file to be opened.
 * @param offset The offset of the read quality.

 * @return Pointer to the new parser object with these filename and
 * offset.
 */
inline Parser* SelectParser(const std::string& filename,
                     int offset = SingleRead::PHRED_OFFSET) {
	return new FastaFastqGzParser(filename, offset);
//	std::string ext = GetExtension(filename);
//	if ((ext == "fastq") || (ext == "fastq.gz") ||
//			(ext == "fasta") || (ext == "fasta.gz")) {
//		return new FastaFastqGzParser(filename, offset);
//	}
//	return NULL;
}

}

#endif /* COMMON_IO_PARSER_HPP */
