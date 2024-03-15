//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once
#include "read_stream.hpp"

#include "threadpool/threadpool.hpp"

namespace io {

template <typename ReadType>
class AsyncReadStream {
    static constexpr size_t BUF_SIZE = 100000;
  public:
    AsyncReadStream(ReadStream<ReadType> stream, ThreadPool::ThreadPool &pool)
            : stream_(std::move(stream)), pool_(pool) {
        init();
    }

    AsyncReadStream(AsyncReadStream&& that) = default;
    ~AsyncReadStream() = default;

    bool eof() { return eof_; }

    AsyncReadStream<ReadType> &operator>>(ReadType &t) {
        VERIFY(!eof());

        auto wait_writing_buffer = [this]() {
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
        };

        if (start_) {
            dispatch_write_job();
            start_ = false;
            wait_writing_buffer();
        }

        t = std::move(read_buffer_[read_pos_++]);

        if (read_pos_ == read_buffer_.size()) {
            wait_writing_buffer();
        }

        return *this;
    }

    void close() {
        if (write_task_.valid())
            write_task_.wait();

        stream_.close();
    }

    bool is_open() {
        return stream_.is_open();
    }

    void reset() {
        if (write_task_.valid())
            write_task_.wait();

        // Reset the stream
        stream_.reset();
        // Restart write jobs
        init();
    }

    constexpr auto && unwrap() {
        return stream_;
    }
    template<class T>
    constexpr auto && recover() {
        return stream_.template recover<T>();
    }

  private:
    void init() {
        read_buffer_.clear();
        write_buffer_.clear();
        read_pos_ = 0;

        read_buffer_.reserve(BUF_SIZE);
        write_buffer_.reserve(BUF_SIZE);
        VERIFY(is_open());

        start_ = true;
        eof_ = stream_.eof();
    }

    void dispatch_write_job() {
        write_task_ =
                pool_.run([this] {
                              while (write_buffer_.size() < BUF_SIZE && !stream_.eof()) {
                                  ReadType r;
                                  stream_ >> r;
                                  write_buffer_.emplace_back(std::move(r));
                              }

                              return !stream_.eof();
                          });
    }

    ReadStream<ReadType> stream_;

    std::vector<ReadType> read_buffer_;
    size_t read_pos_ = 0;
    std::vector<ReadType> write_buffer_;
    bool eof_ = false;
    bool start_ = true;

    std::future<bool> write_task_;
    ThreadPool::ThreadPool &pool_;
};

template<class WrappedStream, typename... Args>
ReadStream<typename WrappedStream::ReadT>
make_async_stream(ThreadPool::ThreadPool &pool, Args&&... args) {
    return AsyncReadStream<typename WrappedStream::ReadT>(WrappedStream(std::forward<Args>(args)...), pool);
}

}
