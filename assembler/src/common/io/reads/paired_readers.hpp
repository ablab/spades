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
         size_t insert_size, bool change_order = false,
         bool use_orientation = true, LibraryOrientation orientation = LibraryOrientation::FR,
         OffsetType offset_type = PhredOffset)
      : insert_size_(insert_size),
        change_order_(change_order),
        use_orientation_(use_orientation),
        changer_(GetOrientationChanger<PairedRead>(orientation)),
        offset_type_(offset_type),
        first_(new FileReadStream(filename1, offset_type_)),
        second_(new FileReadStream(filename2, offset_type_)),
        filename1_(filename1),
        filename2_(filename2){}

  /*
   * Check whether the stream is opened.
   *
   * @return true of the stream is opened and false otherwise.
   */
  /* virtual */ bool is_open() {
    return first_->is_open() && second_->is_open();
  }

  /*
   * Check whether we've reached the end of stream.
   *
   * @return true if the end of stream is reached and false
   * otherwise.
   */
  /* virtual */ bool eof() {

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
  /* virtual */ SeparatePairedReadStream& operator>>(PairedRead& pairedread) {
    SingleRead sr1, sr2;
    (*first_) >> sr1;
    (*second_) >> sr2;

    if (use_orientation_) {
        pairedread = changer_->Perform(PairedRead(sr1, sr2, insert_size_));
    }
    else {
        pairedread = PairedRead(sr1, sr2, insert_size_);
    }

    if (change_order_) {
        pairedread = PairedRead(pairedread.second(), pairedread.first(), insert_size_);
    }

    return *this;
  }

  /*
   * Close the stream.
   */
  /* virtual */ void close() {
    first_->close();
    second_->close();
  }

  /*
   * Close the stream and open it again.
   */
  /* virtual */ void reset() {
    first_->reset();
    second_->reset();
  }

  ReadStreamStat get_stat() const {
    return ReadStreamStat();
  }

 private:

  size_t insert_size_;

  bool change_order_;

  bool use_orientation_;

  std::unique_ptr<OrientationChanger<PairedRead>> changer_;

  /*
   * @variable Quality offset type.
   */
  OffsetType offset_type_;

  /*
   * @variable The first stream (reads from first file).
   */
  std::unique_ptr<ReadStream<SingleRead>> first_;
  /*
   * @variable The second stream (reads from second file).
   */
  std::unique_ptr<ReadStream<SingleRead>> second_;

  //Only for providing information about error for users
  std::string filename1_;
  std::string filename2_;
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
  explicit InterleavingPairedReadStream(const std::string& filename, size_t insert_size, bool change_order = false,
          bool use_orientation = true, LibraryOrientation orientation = LibraryOrientation::FR,
          OffsetType offset_type = PhredOffset)
      : filename_(filename), insert_size_(insert_size),
        change_order_(change_order),
        use_orientation_(use_orientation),
        changer_(GetOrientationChanger<PairedRead>(orientation)),
        offset_type_(offset_type),
        single_(new FileReadStream(filename_, offset_type_)) {}

  /*
   * Check whether the stream is opened.
   *
   * @return true of the stream is opened and false otherwise.
   */
  /* virtual */ bool is_open() {
    return single_->is_open();
  }

  /*
   * Check whether we've reached the end of stream.
   *
   * @return true if the end of stream is reached and false
   * otherwise.
   */
  /* virtual */ bool eof() {
    return single_->eof();
  }

  /*
   * Read PairedRead from stream.
   *
   * @param pairedread The PairedRead that will store read data.
   *
   * @return Reference to this stream.
   */
  /* virtual */ InterleavingPairedReadStream& operator>>(PairedRead& pairedread) {
    SingleRead sr1, sr2;
    (*single_) >> sr1;
    (*single_) >> sr2;

    if (use_orientation_) {
        pairedread = changer_->Perform(PairedRead(sr1, sr2, insert_size_));
    }
    else {
        pairedread = PairedRead(sr1, sr2, insert_size_);
    }

    if (change_order_) {
        pairedread = PairedRead(pairedread.second(), pairedread.first(), insert_size_);
    }

    return *this;
  }

  /*
   * Close the stream.
   */
  /* virtual */ void close() {
    single_->close();
  }

  /*
   * Close the stream and open it again.
   */
  /* virtual */ void reset() {
    single_->reset();
  }

  ReadStreamStat get_stat() const {
        return ReadStreamStat();
  }

 private:
  /*
   * @variable The names of the file which stream reads from.
   */
  std::string filename_;

  size_t insert_size_;

  bool change_order_;

  bool use_orientation_;

  std::unique_ptr<OrientationChanger<PairedRead>> changer_;

  /*
   * @variable Quality offset type.
   */
  OffsetType offset_type_;

  /*
   * @variable The single read stream.
   */
  std::unique_ptr<ReadStream<SingleRead>> single_;

};
}
