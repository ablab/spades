#ifndef __HAMMER_READ_PROCESSOR_HPP__
#define __HAMMER_READ_PROCESSOR_HPP__

#include "read/read.hpp"
#include "read/ireadstream.hpp"

#include "io/mpmc_bounded.hpp"

#include "openmp_wrapper.h"

namespace hammer {
class ReadProcessor {
  unsigned nthreads_;

private:
  template<class Reader, class Op>
  bool RunSingle(Reader &irs, Op &op) {
    while (!irs.eof()) {
      typename Reader::read_type r;
      irs >> r;

      if (op(r))
        return true;
    }

    return false;
  }

  template<class Reader, class Op, class Writer>
  void RunSingle(Reader &irs, Op &op, Writer &writer) {
    while (!irs.eof()) {
      typename Reader::read_type r;
      irs >> r;

      auto res = op(r);

      if (res)
        writer << *res;
    }
  }

public:
  ReadProcessor(unsigned nthreads)
      : nthreads_(nthreads) {}

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

    mpmc_bounded_queue<typename Reader::read_type> in_queue(2*bufsize);

    bool stop = false;
#   pragma omp parallel shared(in_queue, irs, op, stop) num_threads(nthreads_)
    {
#     pragma omp master
      {
        while (!irs.eof()) {
          typename Reader::read_type r;
          irs >> r;

          while (!in_queue.enqueue(r))
            sched_yield();

#         pragma omp flush (stop)
          if (stop)
            break;
        }

        in_queue.close();
      }

      while (1) {
        typename Reader::read_type r;

        if (!in_queue.wait_dequeue(r))
          break;

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

    mpmc_bounded_queue<typename Reader::read_type> in_queue(bufsize), out_queue(bufsize);
#   pragma omp parallel shared(in_queue, out_queue, irs, op, writer) num_threads(nthreads_)
    {
#     pragma omp master
      {
        while (!irs.eof()) {
          typename Reader::read_type r;
          irs >> r;

          // First, try to provide read to the queue. If it's full, never mind.
          bool status = in_queue.enqueue(r);

          // Flush down the output queue
          typename Reader::read_type outr;
          while (out_queue.dequeue(outr))
            writer << outr;

          // If the input queue was originally full, wait until we can insert
          // the read once again.
          if (!status)
            while (!in_queue.enqueue(r))
              sched_yield();
        }

        in_queue.close();
      }

      while (1) {
        typename Reader::read_type r;

        if (!in_queue.wait_dequeue(r))
          break;

        auto res = op(r);
        if (res)
          while (!in_queue.enqueue(*res))
            sched_yield();
      }
    }

    // Flush down the output queue
    typename Reader::read_type outr;
    while (out_queue.dequeue(outr))
      writer << outr;
  }
};

}

#endif // __HAMMER_READ_PROCESSOR_HPP__
