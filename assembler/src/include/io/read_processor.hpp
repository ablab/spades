#ifndef __HAMMER_READ_PROCESSOR_HPP__
#define __HAMMER_READ_PROCESSOR_HPP__

#include "io/mpmc_bounded.hpp"

#include "openmp_wrapper.h"

namespace hammer {
class ReadProcessor {
  static size_t const         cacheline_size = 64;
  typedef char                cacheline_pad_t [cacheline_size];

  unsigned nthreads_;
  cacheline_pad_t pad0;
  size_t read_;
  cacheline_pad_t pad1;
  size_t processed_;
  cacheline_pad_t pad2;

private:
  template<class Reader, class Op>
  bool RunSingle(Reader &irs, Op &op) {
    while (!irs.eof()) {
      typename Reader::ReadT r;
      irs >> r;
      read_ += 1;

      processed_ += 1;
      if (op(r))
        return true;
    }

    return false;
  }

  template<class Reader, class Op, class Writer>
  void RunSingle(Reader &irs, Op &op, Writer &writer) {
    while (!irs.eof()) {
      typename Reader::ReadT r;
      irs >> r;
      read_ += 1;

      auto res = op(r);
      processed_ += 1;

      if (res)
        writer << *res;
    }
  }

public:
  ReadProcessor(unsigned nthreads)
      : nthreads_(nthreads), read_(0), processed_(0) {}

  size_t read() const { return read_; }
  size_t processed() const { return processed_; }

  template<class Reader, class Op>
  bool Run(Reader &irs, Op &op) {
    if (nthreads_ < 2)
      return RunSingle(irs, op);

    // Round nthreads to next power of two
    unsigned bufsize = nthreads_ - 1;
    bufsize = (bufsize >> 1) | bufsize;
    bufsize = (bufsize >> 2) | bufsize;
    bufsize = (bufsize >> 4) | bufsize;
    bufsize = (bufsize >> 8) | bufsize;
    bufsize = (bufsize >> 16) | bufsize;
    bufsize += 1;

    mpmc_bounded_queue<typename Reader::ReadT> in_queue(2*bufsize);

    bool stop = false;
#   pragma omp parallel shared(in_queue, irs, op, stop) num_threads(nthreads_)
    {
#     pragma omp master
      {
        while (!irs.eof()) {
          typename Reader::ReadT r;
          irs >> r;
#         pragma omp atomic
          read_ += 1;

          while (!in_queue.enqueue(r))
            sched_yield();

#         pragma omp flush (stop)
          if (stop)
            break;
        }

        in_queue.close();
      }

      while (1) {
        typename Reader::ReadT r;

        if (!in_queue.wait_dequeue(r))
          break;

#       pragma omp atomic
        processed_ += 1;

        bool res = op(r);
        if (res) {
#         pragma omp atomic
          stop |= res;
        }
      }
    }

#   pragma omp flush(stop)
    return stop;
  }

  template<class Reader, class Op, class Writer>
  void Run(Reader &irs, Op &op, Writer &writer) {
    if (nthreads_ < 2) {
      RunSingle(irs, op, writer);
      return;
    }

    // Round nthreads to next power of two
    unsigned bufsize = nthreads_ - 1;
    bufsize = (bufsize >> 1) | bufsize;
    bufsize = (bufsize >> 2) | bufsize;
    bufsize = (bufsize >> 4) | bufsize;
    bufsize = (bufsize >> 8) | bufsize;
    bufsize = (bufsize >> 16) | bufsize;
    bufsize += 1;

    mpmc_bounded_queue<typename Reader::ReadT> in_queue(bufsize), out_queue(2*bufsize);
#   pragma omp parallel shared(in_queue, out_queue, irs, op, writer) num_threads(nthreads_)
    {
#     pragma omp master
      {
        while (!irs.eof()) {
          typename Reader::ReadT r;
          irs >> r;

          // First, try to provide read to the queue. If it's full, never mind.
          bool status = in_queue.enqueue(r);

          // Flush down the output queue
          typename Reader::ReadT outr;
          while (out_queue.dequeue(outr))
            writer << outr;

          // If the input queue was originally full, wait until we can insert
          // the read once again.
          if (!status)
            while (!in_queue.enqueue(r))
              sched_yield();
        }

        in_queue.close();

        // Flush down the output queue while in master threads.
        typename Reader::ReadT outr;
        while (out_queue.dequeue(outr))
          writer << outr;
      }

      while (1) {
        typename Reader::ReadT r;

        if (!in_queue.wait_dequeue(r))
          break;

        auto res = op(r);
        if (res)
          while (!out_queue.enqueue(*res))
            sched_yield();
      }
    }

    // Flush down the output queue
    typename Reader::ReadT outr;
    while (out_queue.dequeue(outr))
      writer << outr;
  }
};

}

#endif // __HAMMER_READ_PROCESSOR_HPP__
