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
#include "common/io/reader.hpp"

template<typename ReadType>
class MultifileReaderWrapper : public Reader<ReadType> {
 public:
  MultifileReaderWrapper(const vector<typename ReadType::FilenameType>&
                         filenames, size_t distance = 0)
      : filenames_(filenames), distance_(distance) {
    is_open_ = true;
    eof_ = false;
    for (size_t i = 0; i < filenames_.size(); ++i) {
      reader_.push_back(new Reader<ReadType>(filenames_[i], distance));
      is_open_ = is_open_ && reader_[i].is_open();
      eof_ = eof_ || reader_[i].eof();
    }
    current_reader_index_ = 0;
  }

  /*virtual*/~MultifileReaderWrapper() {
    close();
  }

  /*virtual*/MultifileReaderWrapper& operator>>(ReadType& read) {
    if (readers_[current_reader_index_].eof()) {
      ++current_reader_index_;
    }
    if (current_reader_index_ < readers_.size()) {
      (*readers_[current_reader_index_]) >> read;
    } else {
      eof_ = true;
    }
    return (*this);
  }

  /*virtual*/void close() {
    if (is_open_) {
      for (size_t i = 0; i < readers_.size(); ++i) {
        readers_[i].close();
      }
      is_open_ = false;
    }
  }

  /*virtual*/void reset() {
    for (size_t i = 0; i < readers_.size(); ++i) {
      readers_[i].reset();
    }
  }

 private:
  vector<typename ReadType::FilenameType> filenames_;
  vector<Reader<ReadType> > readers_;
  size_t distance_;
  size_t current_reader_index_;
};

#endif /* COMMON_IO_MULTIFILEREADERWRAPPER_HPP_ */
