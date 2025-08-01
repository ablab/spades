//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************
#pragma once

#include "utils/verify.hpp"

#include <boost/noncopyable.hpp>
#include <istream>
#include <memory>
#include <typeinfo>
#include <typeindex>
#include <cstdint>

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

template<class ReadType>
struct ReadStream {
  typedef ReadType ReadT;

  ReadStream() = default;

  template<class T>
  ReadStream(T &&t) noexcept
      : self_(std::make_unique<ReadStreamModel<T>>(std::forward<T>(t))) {}

  bool is_open() const { return self_->is_open(); }
  bool eof() const { return self_->eof(); }
  ReadStream& operator>>(ReadType& read) { (*self_) >> read; return *this; }
  void close() { self_->close(); }
  void reset() { self_->reset(); }

  explicit operator bool() const { return (bool)self_; }

  template<class T>
  constexpr auto &&recover() {
      VERIFY(std::type_index(self_->type_info()) == std::type_index(typeid(T)));

      auto model = static_cast<ReadStreamModel<T>*>(self_.get());
      return model->self_;
  }

 private:
  struct ReadStreamConcept {
    virtual ~ReadStreamConcept() = default;

    /*
     * Check whether the stream is opened.
     * @return true if the stream is opened and false otherwise.
     */
    virtual bool is_open() = 0;

    /*
     * Check whether we've reached the end of stream.
     * @return true if the end of the stream is reached and false
     * otherwise.
     */
    virtual bool eof() = 0;

    /*
     * Read SingleRead or PairedRead from stream (according to ReadType).
     * @param read The SingleRead or PairedRead that will store read data.
     * @return Reference to this stream.
     */
    virtual ReadStreamConcept& operator>>(ReadType& read) = 0;

    /* Close the stream */
    virtual void close() = 0;

    /* Close the stream and open it again */
    virtual void reset() = 0;

    /* Provide the type information for reification */
    virtual const std::type_info &type_info() const = 0;

  };

  template<typename T>
  struct ReadStreamModel : ReadStreamConcept {
    explicit ReadStreamModel(T &&s) noexcept : self_(std::move(s)) {}

    bool is_open() { return self_.is_open(); }
    bool eof() { return self_.eof(); }
    ReadStreamModel& operator>>(ReadType& read) { self_ >> read; return *this; }
    void close() { self_.close(); }
    void reset() { self_.reset(); }
    const std::type_info &type_info() const { return typeid(T); }

    T self_;
  };

  std::unique_ptr<ReadStreamConcept> self_;
};

}
