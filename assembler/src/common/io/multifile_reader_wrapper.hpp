/**
 * @file    multifile_reader_wrapper.hpp
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
 * MultifileReaderWrapper is the class-wrapper that gets data from
 * number of files, given in a constructor.
 */

#ifndef COMMON_IO_MULTIFILEREADERWRAPPER_HPP_
#define COMMON_IO_MULTIFILEREADERWRAPPER_HPP_

#include <vector>
#include "common/io/ireader.hpp"
#include "common/io/reader.hpp"

template<typename ReadType>
class MultifileReaderWrapper : public IReader<ReadType> {
 public:
  /*
   * Default constructor.
   * 
   * @param filenames The names of the files to be opened. Wrapper
   * tries to open all the files of this list and proceeds only those
   * of them which present. The names of non-existing files are
   * ignored. 
   * @param distance Distance between parts of paired reads (or useless
   * parameter when we work with single reads).
   * @param offset The offset of the read quality.
   */
  MultifileReaderWrapper(const vector<typename ReadType::FilenameType>&
                         filenames, size_t distance = 0,
                         int offset = SingleRead::PHRED_OFFSET)
      : filenames_(filenames), distance_(distance), offset_(offset),
        current_reader_index_(0) {
    for (size_t i = 0; i < filenames_.size(); ++i) {
      Reader<ReadType>* reader_ = new Reader<ReadType>(
          filenames_[i], distance_, offset_);
      if (reader_->is_open()) {
        readers_.push_back(reader_);
      } else {
        delete reader_;
      }
    }
  }

  /* 
   * Default destructor.
   */
  /*virtual*/~MultifileReaderWrapper() {
    close();
    for (size_t i = 0; i < readers_.size(); ++i) {
      delete readers_[i];
    }
  }

  /* 
   * Check whether the stream is opened.
   */
  /* virtual */ bool is_open() {
    if (readers_.size() > 0) {
      return readers_[0]->is_open();
    } else {
      return false;
    }
  }

  /* 
   * Check whether we've reached the end of stream.
   */
  /* virtual */ bool eof() {
    if (readers_.size() > 0) {
      return readers_[readers_.size() - 1]->eof();
    } else {
      return true;
    }
  }

  /*
   * Read single or paired read from stream (according to ReadType).
   *
   * @param read The single or paired read that will store read data.
   *
   * @return Reference to this stream.
   */
  /* virtual */ MultifileReaderWrapper& operator>>(ReadType& read) {
    if (readers_.size() > 0) {
      if (readers_[current_reader_index_]->eof()) {
        ++current_reader_index_;
      }
      if (current_reader_index_ < readers_.size()) {
        (*readers_[current_reader_index_]) >> read;
      }
    }
    return (*this);
  }

  /*
   * Close the stream.
   */
  /* virtual */ void close() {
    for (size_t i = 0; i < readers_.size(); ++i) {
      readers_[i]->close();
    }
  }

  /* 
   * Close the stream and open it again.
   */
  /* virtual */ void reset() {
    for (size_t i = 0; i < readers_.size(); ++i) {
      readers_[i]->reset();
    }
  }

 private:
  /*
   * @variable The names of the files which stream read from.
   */
  vector<typename ReadType::FilenameType> filenames_;
  /*
   * @variable Internal stream readers.
   */
  vector<Reader<ReadType>*> readers_;
  /*
   * @variable The distance between two parts of paired read.
   */
  size_t distance_;
  /*
   * @variable Quality offset.
   */
  int offset_;
  /*
   * @variable The index of the file that is currently read from.
   */
  size_t current_reader_index_;
};

#endif /* COMMON_IO_MULTIFILEREADERWRAPPER_HPP_ */
