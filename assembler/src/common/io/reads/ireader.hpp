//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************
//todo rename to reader
#pragma once

#include <boost/noncopyable.hpp>
#include "utils/standard_base.hpp"

namespace io {

struct ReadStreamStat {
    size_t read_count;
    size_t max_len;
    uint64_t total_len;

    ReadStreamStat(): read_count(0), max_len(0), total_len(0) { }

    void write(std::ostream& stream) const {
        stream.write((const char *) &read_count, sizeof(read_count));
        stream.write((const char *) &max_len, sizeof(max_len));
        stream.write((const char *) &total_len, sizeof(total_len));
    }

    void read(std::istream& stream) {
        stream.read((char *) &read_count, sizeof(read_count));
        stream.read((char *) &max_len, sizeof(max_len));
        stream.read((char *) &total_len, sizeof(total_len));
    }

    template<class Read>
    void increase(const Read& read) {
        size_t len = read.size();

        ++read_count;
        if (max_len < len) {
            max_len = len;
        }
        total_len += read.nucl_count();
    }

    void merge(const ReadStreamStat& stat) {
        read_count += stat.read_count;
        if (max_len < stat.max_len) {
            max_len = stat.max_len;
        }
        total_len += stat.total_len;
    }

    bool valid() const {
        return read_count != 0;
    }

};

/**
 * Reader is the interface for all other readers and reader wrappers.
 */
template<typename ReadType>
class ReadStream: boost::noncopyable {
 public:
  typedef ReadType ReadT;

  /*
   * Default destructor.
   */
  virtual ~ReadStream() {}

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
  virtual ReadStream& operator>>(ReadType& read) = 0;

  /*
   * Close the stream.
   */
  virtual void close() = 0;

  /*
   * Close the stream and open it again.
   */
  virtual void reset() = 0;

};

}
