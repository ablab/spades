//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once
#include "io/reads/ireadstream.hpp"

#include "threadpool/threadpool.hpp"

namespace io {

template <typename ReadStreamType>
class AsyncReadStream : public ReadStream<typename ReadStreamType::ReadT> {
    static constexpr size_t BUF_SIZE = 100000;
  public:
    using ReadType = typename ReadStreamType::ReadT;
    
      AsyncReadStream(ReadStream<ReadType> *stream)
              : stream_(stream), pool_(2) {
          init();
      }

      AsyncReadStream(ReadStream<ReadType> &&stream)
              : stream_(stream), pool_(2) {
          init();
      }

      ~AsyncReadStream() { pool_.stop(); }

      bool eof() override { return eof_; }

      AsyncReadStream<ReadStreamType> &operator>>(ReadType &t) override {
          VERIFY(!eof());

          if (read_pos_ < read_buffer_.size()) t = std::move(read_buffer_[read_pos_++]);

          if (read_pos_ == read_buffer_.size()) {
              // Wait for completion of the write task, if any
              bool has_more = false;
              if (write_task_.valid())
                  has_more = write_task_.get();

              // Swap buffers
              std::swap(read_buffer_, write_buffer_);
              read_pos_ = 0;
              write_buffer_.clear();

              // See, if there is still something to read
              eof_ = (read_buffer_.size() == 0);

              // Submit new job
              if (has_more) dispatch_write_job();
          }

          return *this;
    }

    void close() override {
        stream_->close();
        pool_.stop();
    }

    bool is_open() override {
        return stream_->is_open();
    }

    void reset() override {
        // Discard all write jobs
        pool_.stop();
        // Reset the stream
        stream_.reset();
        // Restart write jobs
        init();
    }

  private:
    void init() {
        read_buffer_.clear();
        write_buffer_.clear();
        read_pos_ = 0;
        
        read_buffer_.reserve(BUF_SIZE);
        write_buffer_.reserve(BUF_SIZE);
        VERIFY(is_open());
        dispatch_write_job();
    }

    void dispatch_write_job() {
        write_task_ =
                pool_.run([this] {
                              while (write_buffer_.size() < BUF_SIZE && !stream_->eof()) {
                                  ReadType r;
                                  *stream_ >> r;
                                  write_buffer_.emplace_back(r);
                              }

                              return !stream_->eof();
                          });
    }
    
    std::unique_ptr<ReadStream<ReadType>> stream_;

    std::vector<ReadType> read_buffer_;
    size_t read_pos_ = 0;
    std::vector<ReadType> write_buffer_;
    bool eof_ = false;

    std::future<bool> write_task_;
    ThreadPool::ThreadPool pool_;
};

template<typename ReadStream, typename... Args>
std::unique_ptr<AsyncReadStream<ReadStream>> make_async_stream(Args&&... args) {
    return std::make_unique<AsyncReadStream<ReadStream>>(new ReadStream{std::forward<Args>(args)...});
}


}
