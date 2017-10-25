//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <string>
#include "ireader.hpp"
#include "paired_read.hpp"
#include "file_reader.hpp"
#include "orientation.hpp"

namespace io {

class SeparatePairedReadStream : public ReadStream<PairedRead> {
 public:
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
         OffsetType offset_type = PhredOffset)
      : insert_size_(insert_size),
        first_(new FileReadStream(filename1, offset_type)),
        second_(new FileReadStream(filename2, offset_type)),
        filename1_(filename1),
        filename2_(filename2){}

  /*
   * Check whether the stream is opened.
   *
   * @return true of the stream is opened and false otherwise.
   */
  bool is_open() override {
    return first_->is_open() && second_->is_open();
  }

  /*
   * Check whether we've reached the end of stream.
   *
   * @return true if the end of stream is reached and false
   * otherwise.
   */
  bool eof() override {

    if (first_->eof() != second_->eof()) {
        if (first_->eof()) {
            ERROR("The number of right read-pairs is larger than the number of left read-pairs");
        } else {
            ERROR("The number of left read-pairs is larger than the number of right read-pairs");
        }
        FATAL_ERROR("Unequal number of read-pairs detected in the following files: " << filename1_ << "  " << filename2_ << "");
    }
    return first_->eof();
  }

  /*
   * Read PairedRead from stream.
   *
   * @param pairedread The PairedRead that will store read data.
   *
   * @return Reference to this stream.
   */
  SeparatePairedReadStream& operator>>(PairedRead& pairedread) override {
    pairedread.set_orig_insert_size(insert_size_);
    (*first_) >> pairedread.first();
    (*second_) >> pairedread.second();
    return *this;
  }

  /*
   * Close the stream.
   */
  void close() override {
    first_->close();
    second_->close();
  }

  /*
   * Close the stream and open it again.
   */
  void reset() override {
    first_->reset();
    second_->reset();
  }

 private:
  const size_t insert_size_;

  /*
   * @variable The first stream (reads from first file).
   */
  const std::unique_ptr<ReadStream<SingleRead>> first_;
  /*
   * @variable The second stream (reads from second file).
   */
  const std::unique_ptr<ReadStream<SingleRead>> second_;

  //Only for providing information about error for users
  const std::string filename1_;
  const std::string filename2_;
};

class InterleavingPairedReadStream : public ReadStream<PairedRead> {
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
                                        OffsetType offset_type = PhredOffset)
      : filename_(filename), insert_size_(insert_size),
        single_(new FileReadStream(filename_, offset_type)) {}

  /*
   * Check whether the stream is opened.
   *
   * @return true of the stream is opened and false otherwise.
   */
  bool is_open() override {
    return single_->is_open();
  }

  /*
   * Check whether we've reached the end of stream.
   *
   * @return true if the end of stream is reached and false
   * otherwise.
   */
  bool eof() override {
    return single_->eof();
  }

  /*
   * Read PairedRead from stream.
   *
   * @param pairedread The PairedRead that will store read data.
   *
   * @return Reference to this stream.
   */
  InterleavingPairedReadStream& operator>>(PairedRead& pairedread) override {
    pairedread.set_orig_insert_size(insert_size_);
    (*single_) >> pairedread.first();
    (*single_) >> pairedread.second();
    return *this;
  }

  /*
   * Close the stream.
   */
  void close() override {
    single_->close();
  }

  /*
   * Close the stream and open it again.
   */
  void reset() override {
    single_->reset();
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
  const std::unique_ptr<ReadStream<SingleRead>> single_;

};
}
