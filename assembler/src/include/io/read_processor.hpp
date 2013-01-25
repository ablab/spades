#ifndef __HAMMER_READ_PROCESSOR_HPP__
#define __HAMMER_READ_PROCESSOR_HPP__

#include "read/read.hpp"
#include "read/ireadstream.hpp"

#include <gcl/buffer_queue.h>

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
    gcl::buffer_queue<typename Reader::read_type> in_queue(nthreads_);

    if (nthreads_ < 2)
      return RunSingle(irs, op);

    bool stop = false;
#   pragma omp parallel shared(in_queue, irs, op, stop) num_threads(nthreads_)
    {
#     pragma omp master
      {
        while (!irs.eof()) {
          typename Reader::read_type r;
          irs >> r;

          auto status = in_queue.wait_push(r);
          VERIFY(status == gcl::queue_op_status::success);

#         pragma omp flush(stop)
          if (stop)
            break;
        }

        in_queue.close();
      }

      while (!in_queue.is_closed()) {
      typename Reader::read_type r;
        if (in_queue.wait_pop(r) == gcl::queue_op_status::closed)
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
    gcl::buffer_queue<typename Reader::read_type> in_queue(nthreads_), out_queue(nthreads_);

    if (nthreads_ < 2) {
      RunSingle(irs, op, writer);
      return;
    }

#   pragma omp parallel shared(in_queue, out_queue, irs, op, writer) num_threads(nthreads_)
    {
#     pragma omp master
      {
        while (!irs.eof()) {
          typename Reader::read_type r;
          irs >> r;

          // First, try to provide read to the queue. If it's full, never mind.
          auto status = in_queue.try_push(r);
          VERIFY(status == gcl::queue_op_status::success ||
                 status == gcl::queue_op_status::full);

          // Flush down the output queue
          while (!out_queue.is_empty())
            writer << out_queue.value_pop();

          // If the input queue was originally full, wait until we can insert
          // the read once again.
          if (status == gcl::queue_op_status::full)
            in_queue.push(r);
        }

        in_queue.close();
      }

      while (!in_queue.is_closed()) {
        typename Reader::read_type r;
        if (in_queue.wait_pop(r) == gcl::queue_op_status::closed)
          break;

        auto res = op(r);
        if (res)
          out_queue.push(*res);
      }
    }

    // Flush down the output queue
    while (!out_queue.is_empty())
      writer << out_queue.value_pop();
  }

};

}

#endif // __HAMMER_READ_PROCESSOR_HPP__
