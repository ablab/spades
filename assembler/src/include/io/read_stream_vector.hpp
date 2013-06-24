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

#pragma once

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

//todo check destroy_readers logic and usages
template <class Reader>
class ReadStreamVector : boost::noncopyable {

 private:
  mutable std::vector<Reader*> streams_;
  bool destroy_readers_;

  ReadStreamVector& operator=(const ReadStreamVector& that);

 public:

  typedef Reader ReaderType;

  ReadStreamVector(const std::vector<Reader*>& streams, bool destroy_readers = true): streams_(streams.size()), destroy_readers_(destroy_readers) {
    std::copy(streams.begin(), streams.end(), streams_.begin());
  }

  ReadStreamVector(bool destroy_readers = true): destroy_readers_(destroy_readers) {
  }

  ReadStreamVector(Reader* stream): streams_(1, stream), destroy_readers_(true) {
  }

  ReadStreamVector(Reader& stream): streams_(1, &stream), destroy_readers_(false) {
  }

  std::vector<Reader*>& get() {
      destroy_readers_ = false;
      return streams_;
  }

//  ReadStreamVector(const ReadStreamVector& that): streams_(that.streams_), destroy_readers_(false) {
//  }

  ~ReadStreamVector() {
    if (destroy_readers_) {
      for (size_t i = 0; i < streams_.size(); ++i) {
        delete streams_[i];
      }
    }
  }

  class iterator: public std::iterator<std::input_iterator_tag, Reader> {
    typedef typename std::vector<Reader*>::iterator vec_it;
    vec_it it_;
   public:

    iterator(vec_it it) : it_(it) {
    }

    void operator++ () {
        ++it_;
    }

    bool operator== (const iterator& that) {
        return it_ == that.it_;
    }

    bool operator!= (const iterator& that) {
        return it_ != that.it_;
    }

    Reader& operator*() {
        return *(*it_);
    }
  };

  class const_iterator: public std::iterator<std::input_iterator_tag, Reader> {
    typedef typename std::vector<Reader*>::iterator vec_it;
    vec_it it_;
   public:

    const_iterator(vec_it it) : it_(it) {
    }

    void operator++ () {
        ++it_;
    }

    bool operator== (const const_iterator& that) {
        return it_ == that.it_;
    }

    bool operator!= (const const_iterator& that) {
        return it_ != that.it_;
    }

    Reader& operator*() {
        return *(*it_);
    }
  };

  Reader& operator[](size_t i) {
    return *streams_.at(i);
  }

  Reader& back() {
    return *streams_.back();
  }

  size_t size() const {
    return streams_.size();
  }

  bool eof() const {
    return ParallelStreamEOF(streams_);
  }

  iterator begin() {
    return iterator(streams_.begin());
  }

  iterator end() {
    return iterator(streams_.end());
  }

  const_iterator begin() const {
    return iterator(streams_.begin());
  }

  const_iterator end() const {
    return iterator(streams_.end());
  }

  void push_back(Reader* stream) {
      streams_.push_back(stream);
  }

  void reset() {
    for (auto I = streams_.begin(), E = streams_.end(); I != E; ++I)
      (*I)->reset();
  }

  void release() {
      destroy_readers_ = false;
  }

  const std::vector< Reader * >& get() const {
      return streams_;
  }

};

}
