//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*  Multi-consumer/multi-producer bounded queue

Copyright (c) 2011 Dmitry Vyukov. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list
of conditions and the following disclaimer in the documentation and/or other materials
provided with the distribution.

THIS SOFTWARE IS PROVIDED BY DMITRY VYUKOV "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL DMITRY VYUKOV OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include <ciso646>

#if __GNUC__ > 4 || (__GNUC__ >= 4 && __GNUC_MINOR__ >= 5) || _LIBCPP_VERSION

#include <atomic>

#else
#include <cstdatomic>
#endif

#include <cstring>
#include <unistd.h>

template<typename T>
class mpmc_bounded_queue {
public:
    mpmc_bounded_queue(size_t buffer_size)
            : buffer_(new cell_t[buffer_size]), buffer_mask_(buffer_size - 1) {
        assert((buffer_size >= 2) && ((buffer_size & (buffer_size - 1)) == 0));
        for (size_t i = 0; i != buffer_size; i += 1)
            buffer_[i].sequence_.store(i, std::memory_order_relaxed);
        enqueue_pos_.store(0, std::memory_order_relaxed);
        dequeue_pos_.store(0, std::memory_order_relaxed);
        closed_.store(false, std::memory_order_relaxed);
    }

    ~mpmc_bounded_queue() {
        delete[] buffer_;
    }

    bool is_closed() const {
        return closed_.load(std::memory_order_relaxed);
    }

    void close() {
        closed_.store(true, std::memory_order_release);
    }

    bool enqueue(T const &data) {
        if (is_closed())
            return false;

        cell_t *cell;
        size_t pos = enqueue_pos_.load(std::memory_order_relaxed);
        for (; ;) {
            cell = &buffer_[pos & buffer_mask_];
            size_t seq = cell->sequence_.load(std::memory_order_acquire);
            intptr_t dif = (intptr_t) seq - (intptr_t) pos;
            if (dif == 0) {
                if (enqueue_pos_.compare_exchange_weak(pos, pos + 1, std::memory_order_relaxed))
                    break;
            } else if (dif < 0)
                return false;
            else
                pos = enqueue_pos_.load(std::memory_order_relaxed);
        }

        cell->data_ = data;
        cell->sequence_.store(pos + 1, std::memory_order_release);

        return true;
    }

    bool enqueue(T &&data) {
        if (is_closed())
            return false;

        cell_t *cell;
        size_t pos = enqueue_pos_.load(std::memory_order_relaxed);
        for (; ;) {
            cell = &buffer_[pos & buffer_mask_];
            size_t seq = cell->sequence_.load(std::memory_order_acquire);
            intptr_t dif = (intptr_t) seq - (intptr_t) pos;
            if (dif == 0) {
                if (enqueue_pos_.compare_exchange_weak(pos, pos + 1, std::memory_order_relaxed))
                    break;
            } else if (dif < 0)
                return false;
            else
                pos = enqueue_pos_.load(std::memory_order_relaxed);
        }

        cell->data_ = std::move(data);
        cell->sequence_.store(pos + 1, std::memory_order_release);

        return true;
    }

    bool dequeue(T &data) {
        cell_t *cell;
        size_t pos = dequeue_pos_.load(std::memory_order_relaxed);
        for (; ;) {
            cell = &buffer_[pos & buffer_mask_];
            size_t seq = cell->sequence_.load(std::memory_order_acquire);
            intptr_t dif = (intptr_t) seq - (intptr_t) (pos + 1);
            if (dif == 0) {
                if (dequeue_pos_.compare_exchange_weak(pos, pos + 1, std::memory_order_relaxed))
                    break;
            } else if (dif < 0)
                return false;
            else
                pos = dequeue_pos_.load(std::memory_order_relaxed);
        }

        data = std::move(cell->data_);
        cell->sequence_.store(pos + buffer_mask_ + 1, std::memory_order_release);

        return true;
    }

    bool wait_dequeue(T &data) {
        bool res = false;
        do {
            res = dequeue(data);
            if (!res)
                usleep(1);
        } while (!res && !is_closed());

        return res;
    }

private:
    struct cell_t {
        std::atomic<size_t> sequence_;
        T data_;
    };

    static size_t const cacheline_size = 64;
    typedef char cacheline_pad_t[cacheline_size];

    cacheline_pad_t pad0_;
    cell_t *const buffer_;
    size_t const buffer_mask_;
    cacheline_pad_t pad1_;
    std::atomic<size_t> enqueue_pos_;
    cacheline_pad_t pad2_;
    std::atomic<size_t> dequeue_pos_;
    cacheline_pad_t pad3_;
    std::atomic<bool> closed_;
    cacheline_pad_t pad4_;

    mpmc_bounded_queue(mpmc_bounded_queue const &);

    void operator=(mpmc_bounded_queue const &);
};
