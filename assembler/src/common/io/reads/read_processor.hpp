//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef __HAMMER_READ_PROCESSOR_HPP__
#define __HAMMER_READ_PROCESSOR_HPP__

#include "io/reads/mpmc_bounded.hpp"

#include "utils/parallel/openmp_wrapper.h"

#pragma GCC diagnostic push
#ifdef __clang__
#pragma clang diagnostic ignored "-Wunused-private-field"
#endif
namespace hammer {
class ReadProcessor {
    static size_t constexpr cacheline_size = 64;
    typedef char cacheline_pad_t[cacheline_size];

    unsigned nthreads_;
    cacheline_pad_t pad0;
    size_t read_;
    cacheline_pad_t pad1;
    size_t processed_;
    cacheline_pad_t pad2;

private:
    template<class Reader, class Op>
    bool RunSingle(Reader &irs, Op &op) {
        using ReadPtr = std::unique_ptr<typename Reader::ReadT>;

        while (!irs.eof()) {
            ReadPtr r = ReadPtr(new typename Reader::ReadT) ;
            irs >> *r;
            read_ += 1;

            processed_ += 1;
            if (op(std::move(r))) // Pass ownership of read down to processor
                return true;
        }

        return false;
    }

    template<class Reader, class Op, class Writer>
    void RunSingle(Reader &irs, Op &op, Writer &writer) {
        using ReadPtr = std::unique_ptr<typename Reader::ReadT>;

        while (!irs.eof()) {
            ReadPtr r = ReadPtr(new typename Reader::ReadT) ;
            irs >> *r;
            read_ += 1;

            auto res = op(std::move(r)); // Pass ownership of read down to processor
            processed_ += 1;

            if (res)
                writer << *res;
        }
    }

public:
    ReadProcessor(unsigned nthreads)
            : nthreads_(nthreads), read_(0), processed_(0) { }

    size_t read() const { return read_; }

    size_t processed() const { return processed_; }

    template<class Reader, class Op>
    bool Run(Reader &irs, Op &op) {
        using ReadPtr = std::unique_ptr<typename Reader::ReadT>;

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

        mpmc_bounded_queue<ReadPtr> in_queue(2 * bufsize);

        bool stop = false;
#   pragma omp parallel shared(in_queue, irs, op, stop) num_threads(nthreads_)
        {
#     pragma omp master
            {
                while (!irs.eof()) {
                    ReadPtr r = ReadPtr(new typename Reader::ReadT) ;
                    irs >> *r;
#         pragma omp atomic
                    read_ += 1;

                    while (!in_queue.enqueue(std::move(r)))
                        sched_yield();

#         pragma omp flush (stop)
                    if (stop)
                        break;
                }

                in_queue.close();
            }

            while (1) {
                ReadPtr r;

                if (!in_queue.wait_dequeue(r))
                    break;

#       pragma omp atomic
                processed_ += 1;

                bool res = op(std::move(r));
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
        using ReadPtr = std::unique_ptr<typename Reader::ReadT>;

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

        mpmc_bounded_queue<ReadPtr> in_queue(bufsize), out_queue(2 * bufsize);
#   pragma omp parallel shared(in_queue, out_queue, irs, op, writer) num_threads(nthreads_)
        {
#     pragma omp master
            {
                while (!irs.eof()) {
                    ReadPtr r = ReadPtr(new typename Reader::ReadT) ;
                    irs >> *r;

                    // First, try to provide read to the queue. If it's full, never mind.
                    bool status = in_queue.enqueue(std::move(r));

                    // Flush down the output queue
                    ReadPtr outr;
                    while (out_queue.dequeue(outr))
                        writer << *outr;

                    // If the input queue was originally full, wait until we can insert
                    // the read once again.
                    if (!status)
                        while (!in_queue.enqueue(std::move(r)))
                            sched_yield();
                }

                in_queue.close();

                // Flush down the output queue while in master threads.
                ReadPtr outr;
                while (out_queue.dequeue(outr))
                    writer << *outr;
            }

            while (1) {
                ReadPtr r;

                if (!in_queue.wait_dequeue(r))
                    break;

                auto res = op(std::move(r));
                if (res)
                    while (!out_queue.enqueue(std::move(res)))
                        sched_yield();
            }
        }

        // Flush down the output queue
        ReadPtr outr;
        while (out_queue.dequeue(outr))
            writer << *outr;
    }
};

#pragma GCC diagnostic pop

}

#endif // __HAMMER_READ_PROCESSOR_HPP__
