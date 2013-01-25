//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/**
 * @file    ireader.hpp
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
 * IReader is the interface for all other readers and reader wrappers.
 */

#ifndef COMMON_IO_IREADER_HPP_
#define COMMON_IO_IREADER_HPP_

#include <boost/noncopyable.hpp>

namespace io {


struct ReadStat {

    size_t read_count_;

    size_t max_len_;

    u_int64_t total_len_;


    ReadStat(): read_count_(0), max_len_(0), total_len_(0) {
    }

    void write(std::ostream& stream) const {
        stream.write((const char *) &read_count_, sizeof(read_count_));
        stream.write((const char *) &max_len_, sizeof(max_len_));
        stream.write((const char *) &total_len_, sizeof(total_len_));
    }

    void read(std::istream& stream) {
        stream.read((char *) &read_count_, sizeof(read_count_));
        stream.read((char *) &max_len_, sizeof(max_len_));
        stream.read((char *) &total_len_, sizeof(total_len_));
    }

    template<class Read>
    void increase(const Read& read) {
        size_t len = read.size();

        ++read_count_;
        if (max_len_ < len) {
            max_len_ = len;
        }
        total_len_ += read.nucl_count();
    }

    void merge(const ReadStat& stat) {
        read_count_ += stat.read_count_;
        if (max_len_ < stat.max_len_) {
            max_len_ = stat.max_len_;
        }
        total_len_ += stat.total_len_;
    }

    bool valid() const {
        return read_count_ != 0;
    }

};



template<typename ReadType>
class IReader: boost::noncopyable {
 public:


  typedef ReadType read_type;


  /* 
   * Default destructor.
   */
  virtual ~IReader() {}

  /* 
   * Check whether the stream is opened.
   *
   * @return true if the stream is opened and false otherwise.
   */
  virtual bool is_open() = 0;

  /* 
   * Check whether we've reached the end of stream.
   *
   * @return true if the end of the stream is reached and false
   * otherwise.
   */
  virtual bool eof() = 0;

  /*
   * Read SingleRead or PairedRead from stream (according to ReadType).
   *
   * @param read The SingleRead or PairedRead that will store read data.
   *
   * @return Reference to this stream.
   */
  virtual IReader& operator>>(ReadType& read) = 0;

  /*
   * Close the stream.
   */
  virtual void close() = 0;

  /* 
   * Close the stream and open it again.
   */
  virtual void reset() = 0;


  virtual ReadStat get_stat() const = 0;

};

}

#endif /* COMMON_IO_IREADER_HPP_ */
