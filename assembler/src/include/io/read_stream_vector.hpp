//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * read_stream_vector.hpp
 *
 *  Created on: Jul 24, 2012
 *      Author: andrey
 */

#ifndef READ_STREAM_VECTOR_HPP_
#define READ_STREAM_VECTOR_HPP_


namespace io {


template <class Reader>
bool ParallelStreamEOF(const std::vector< Reader* >& streams) {
  for (size_t i = 0; i < streams.size(); ++i) {
    if (!streams[i]->eof()) {
      return false;
    }
  }
  return true;
}


template <class Reader>
class ReadStreamVector {

 private:
  std::vector< Reader * > streams_;

  bool destroy_readers_;

  ReadStreamVector& operator=(const ReadStreamVector& that);

 public:

  ReadStreamVector(std::vector< Reader *> streams, bool destroy_readers = true): streams_(streams.size()), destroy_readers_(destroy_readers) {
    for (size_t i = 0; i < streams.size(); ++i) {
      streams_[i] = streams[i];
    }
  }

  ReadStreamVector(Reader* stream): streams_(1, stream), destroy_readers_(false) {
  }

  ReadStreamVector(const ReadStreamVector& that): streams_(that.streams_), destroy_readers_(false) {
  }

  ~ReadStreamVector() {
    if (destroy_readers_) {
      for (size_t i = 0; i < streams_.size(); ++i) {
        delete streams_[i];
      }
    }
  }

  Reader& operator[](size_t i) const {
    return *streams_.at(i);
  }

  Reader& back() const {
    return *streams_.back();
  }

  size_t size() const {
    return streams_.size();
  }

  bool eof() const {
    return ParallelStreamEOF(streams_);
  }

  void reset() const {
    for (auto I = streams_.begin(), E = streams_.end(); I != E; ++I)
      (*I)->reset();
  }

};

}

#endif /* READ_STREAM_VECTOR_HPP_ */
