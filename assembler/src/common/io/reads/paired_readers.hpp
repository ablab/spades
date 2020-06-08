//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "read_stream.hpp"
#include "paired_read.hpp"
#include "orientation.hpp"

#include <string>

namespace ThreadPool {
class ThreadPool;
};

namespace io {

class SeparatePairedReadStream {
 public:
  typedef PairedRead ReadT;
   
  /*
   * Default constructor.
   *
   * @param filename The pair that contains the names of two files to
   * be opened.
   * @param distance Distance between parts of PairedReads.
   * @param offset The offset of the read quality.
   */
  explicit SeparatePairedReadStream(const std::string& filename1, const std::string& filename2,
                                    size_t insert_size,
                                    FileReadFlags flags = FileReadFlags(),
                                    ThreadPool::ThreadPool *pool = nullptr);

  /*
   * Check whether the stream is opened.
   *
   * @return true of the stream is opened and false otherwise.
   */
  bool is_open() {
    return first_.is_open() && second_.is_open();
  }

  /*
   * Check whether we've reached the end of stream.
   *
   * @return true if the end of stream is reached and false
   * otherwise.
   */
  bool eof();

  /*
   * Read PairedRead from stream.
   *
   * @param pairedread The PairedRead that will store read data.
   *
   * @return Reference to this stream.
   */
  SeparatePairedReadStream& operator>>(PairedRead& pairedread);

  /*
   * Close the stream.
   */
  void close() {
    first_.close();
    second_.close();
  }

  /*
   * Close the stream and open it again.
   */
  void reset() {
    first_.reset();
    second_.reset();
  }

 private:
  const size_t insert_size_;

  /*
   * @variable The first stream (reads from first file).
   */
  ReadStream<SingleRead> first_;
  /*
   * @variable The second stream (reads from second file).
   */
  ReadStream<SingleRead> second_;

  //Only for providing information about error for users
  const std::string filename1_;
  const std::string filename2_;
};

class InterleavingPairedReadStream {
 public:
  /*
   * Default constructor.
   *
   * @param filename Single file
   * @param distance Distance between parts of PairedReads.
   * @param offset The offset of the read quality.
   */
  explicit InterleavingPairedReadStream(const std::string& filename,
                                        size_t insert_size,
                                        FileReadFlags flags = FileReadFlags(),
                                        ThreadPool::ThreadPool *pool = nullptr);
  /*
   * Check whether the stream is opened.
   *
   * @return true of the stream is opened and false otherwise.
   */
  bool is_open() {
    return single_.is_open();
  }

  /*
   * Check whether we've reached the end of stream.
   *
   * @return true if the end of stream is reached and false
   * otherwise.
   */
  bool eof() {
    return single_.eof();
  }

  /*
   * Read PairedRead from stream.
   *
   * @param pairedread The PairedRead that will store read data.
   *
   * @return Reference to this stream.
   */
  InterleavingPairedReadStream& operator>>(PairedRead& pairedread);

  /*
   * Close the stream.
   */
  void close() {
    single_.close();
  }

  /*
   * Close the stream and open it again.
   */
  void reset() {
    single_.reset();
  }

 private:
  /*
   * @variable The names of the file which stream reads from.
   */
  const std::string filename_;
  const size_t insert_size_;

  /*
   * @variable The single read stream.
   */
  ReadStream<SingleRead> single_;
};
}
