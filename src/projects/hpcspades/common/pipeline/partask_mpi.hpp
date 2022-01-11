//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2019-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io/binary/binary.hpp"

#include "utils/parallel/openmp_wrapper.h"
#include "utils/verify.hpp"
#include "utils/logger/logger.hpp"
#include "utils/stl_utils.hpp"
#include "utils/perf/timetracer.hpp"

#include <algorithm>
#include <functional>
#include <limits>
#include <memory>
#include <vector>
#include <numeric>
#include <streambuf>
#include <tuple>

#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <mpi.h>

#define ASSERT_MAIN_THREAD VERIFY(omp_get_thread_num() == 0)

namespace partask {

// Sasha's own idea, motivated by std::as_const function
template <class T>
constexpr std::add_const_t<T> *as_const(T *p) noexcept {
    return p;
}

// Motivated by as_const
template <typename T>
constexpr T &as_non_const(const T &t) noexcept {
    return const_cast<T &>(t);
}

template <typename T>
constexpr T *as_non_const(const T *p) noexcept {
    return const_cast<T *>(p);
}

template <typename T>
auto wrap_into_ptr(T&& t) {
    return new std::remove_reference_t<T>(std::forward<T>(t));
}

template <typename T>
auto wrap_into_shared(T&& t) {
    return std::shared_ptr<std::remove_reference_t<T>>(wrap_into_ptr(std::forward<T>(t)));
}

template <typename T>
auto wrap_into_unique(T&& t) {
    return std::unique_ptr<std::remove_reference_t<T>>(wrap_into_ptr(std::forward<T>(t)));
}


template <typename T>
std::unique_ptr<T> as_unique(T *p) {
    return std::unique_ptr<T>(p);
}

template <typename T>
std::shared_ptr<T> as_shared(T *p) {
    return std::shared_ptr<T>(p);
}

// Motivated by declval, the same stuff but returns lvalue reference
template <class T>
typename std::add_lvalue_reference<T>::type declref() noexcept;

inline int world_size() {
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    return world_size;
}

inline int world_rank() {
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    return world_rank;
}

inline bool master(int rank = 0) {
    return world_rank() == rank;
}

inline bool worker(int rank = 0) {
    return world_rank() != rank;
}

inline bool initialized() {
    int flag;
    MPI_Initialized(&flag);
    return flag;
}

inline void barrier() {
    ASSERT_MAIN_THREAD;
    TIME_TRACE_SCOPE("partask::barrier");
    static size_t count = 0;
    DEBUG("barrier() called " << count << " times");
    ++count;
    int ret = MPI_Barrier(MPI_COMM_WORLD);
    VERIFY(ret == MPI_SUCCESS);
}

const size_t MPI_MAX_COUNT = 1 << 30; // Should be <= MAX_INT

inline void membroadcast(void *p, size_t count, int root = 0) {
    ASSERT_MAIN_THREAD;
    TIME_TRACE_SCOPE("partask::membroadcast");
    static size_t call_count = 0;
    DEBUG("membroadcast() called " << call_count << " times");
    ++call_count;

    char *cp = reinterpret_cast<char*>(p);
    while (count) {
        size_t block_size = std::min(count, MPI_MAX_COUNT);
        int ret = MPI_Bcast(cp, static_cast<int>(block_size), MPI_BYTE, root, MPI_COMM_WORLD);
        VERIFY(ret == MPI_SUCCESS);
        cp += block_size;
        count -= block_size;
    }
}

inline void memsend(const void *p, size_t count, int rank, int tag = 0) {
    ASSERT_MAIN_THREAD;
    TIME_TRACE_SCOPE("partask::memsend");
    char *cp = reinterpret_cast<char*>(const_cast<void*>(p));
    while (count) {
        size_t block_size = std::min(count, MPI_MAX_COUNT);
        int ret = MPI_Send(cp, static_cast<int>(block_size), MPI_BYTE, rank, tag, MPI_COMM_WORLD);
        VERIFY(ret == MPI_SUCCESS);
        cp += block_size;
        count -= block_size;
    }
}

inline void memrecv(void *p, size_t count, int rank, int tag = MPI_ANY_TAG) {
    ASSERT_MAIN_THREAD;
    TIME_TRACE_SCOPE("partask::memrecv");
    char *cp = reinterpret_cast<char*>(p);
    while (count) {
        size_t block_size = std::min(count, MPI_MAX_COUNT);
        int ret = MPI_Recv(cp, static_cast<int>(block_size), MPI_BYTE, rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        VERIFY(ret == MPI_SUCCESS);
        cp += block_size;
        count -= block_size;
    }
}

template <typename T>
std::enable_if_t<std::is_pod<T>::value && !std::is_pointer<T>::value> broadcast(T &v, int root = 0) {
    membroadcast(&v, sizeof(T), root);
}

template <typename T, size_t N>
std::enable_if_t<std::is_pod<T>::value && !std::is_pointer<T>::value> broadcast(std::array<T, N> &v, int root = 0) {
    membroadcast(v.deta(), N * sizeof(T), root);
}

template <typename T>
std::enable_if_t<std::is_pod<T>::value && !std::is_pointer<T>::value> send(const T &v, int rank, int tag) {
    memsend(&v, sizeof(T), rank, tag);
}

template <typename T>
std::enable_if_t<std::is_pod<T>::value && !std::is_pointer<T>::value> recv(T &v, int rank, int tag = MPI_ANY_TAG) {
    memrecv(&v, sizeof(T), rank, tag);
}

template <typename T, size_t N>
std::enable_if_t<std::is_pod<T>::value && !std::is_pointer<T>::value> send(const std::array<T, N> &v, int rank, int tag = 0) {
    memsend(v.data(), N * sizeof(T), rank, tag);
}

template <typename T, size_t N>
std::enable_if_t<std::is_pod<T>::value && !std::is_pointer<T>::value> recv(std::array<T, N> &v, int rank, int tag = MPI_ANY_TAG) {
    memrecv(v.data(), N * sizeof(T), rank, tag);
}

// TODO Add send/recv/broadcast overloadings for vectors of PODs

template <typename T>
inline MPI_Datatype mpi_datatype();

// We use MPI_BYTE for char since it's more safe
// Using MPI_BYTE prevents MPI from performing any kind of representation conversion
template <>
inline MPI_Datatype mpi_datatype<void>() {
    return MPI_BYTE;
}

template <>
inline MPI_Datatype mpi_datatype<char>() {
    return MPI_CHAR;
}

template <>
inline MPI_Datatype mpi_datatype<signed char>() {
    return MPI_SIGNED_CHAR;
}

template <>
inline MPI_Datatype mpi_datatype<unsigned char>() {
    return MPI_UNSIGNED_CHAR;
}

template <>
inline MPI_Datatype mpi_datatype<short>() {
    return MPI_SHORT;
}

template <>
inline MPI_Datatype mpi_datatype<int>() {
    return MPI_INT;
}

template <>
inline MPI_Datatype mpi_datatype<long>() {
    return MPI_LONG;
}

template <>
inline MPI_Datatype mpi_datatype<long long>() {
    return MPI_LONG_LONG;
}

template <>
inline MPI_Datatype mpi_datatype<float>() {
    return MPI_FLOAT;
}

template <>
inline MPI_Datatype mpi_datatype<double>() {
    return MPI_DOUBLE;
}

template <>
inline MPI_Datatype mpi_datatype<long double>() {
    return MPI_LONG_DOUBLE;
}

template <>
inline MPI_Datatype mpi_datatype<unsigned short>() {
    return MPI_UNSIGNED_SHORT;
}

template <>
inline MPI_Datatype mpi_datatype<unsigned>() {
    return MPI_UNSIGNED;
}

template <>
inline MPI_Datatype mpi_datatype<unsigned long>() {
    return MPI_UNSIGNED_LONG;
}

template <>
inline MPI_Datatype mpi_datatype<unsigned long long>() {
    return MPI_UNSIGNED_LONG_LONG;
}

template <typename T>
void allreduce(T *recvbuf, size_t count, MPI_Op op) {
    ASSERT_MAIN_THREAD;
    TIME_TRACE_SCOPE("partask::allreduce");
    DEBUG("allreduce started for " << count << " objects of type " << typeid(T).name());
    using NoneVoidT = std::conditional_t<std::is_void<T>::value, char, T>;
    NoneVoidT *crecvbuf = reinterpret_cast<NoneVoidT*>(recvbuf);
    while (count) {
        size_t block_size = std::min(count, MPI_MAX_COUNT);
        int ret = MPI_Allreduce(MPI_IN_PLACE, crecvbuf, static_cast<int>(block_size), mpi_datatype<T>(), op, MPI_COMM_WORLD);
        VERIFY(ret == MPI_SUCCESS);
        crecvbuf += block_size;
        count -= block_size;
    }
    DEBUG("allreduce finished");
}

namespace detail {

template <typename T, bool IsReference>
struct wrap {};

template <typename T>
struct wrap<T, false> {
    const char *name() const { return typeid(T).name(); }

    T get() { return std::move(t_); }

    template <typename TT>
    wrap(TT &&t) : t_{std::forward<TT>(t)} {}

    T t_;
};

template <typename T>
struct wrap<T, true> {
    const char *name() const { return typeid(T).name(); }

    T &get() { return t_; }

    wrap(T &t) : t_{t} {}

    T t_;
};

template <>
struct wrap<void, false> {
    const char *name() const { return "VOID"; }

    void get() const { return; }
};

using wrap_void = wrap<void, false>;

template <typename T>
auto operator,(T &&t, wrap<void, false>) {
    static_assert(!std::is_reference<T>::value, "TODO Add message");
    DEBUG("T&& called");
    return wrap<T, false>(std::forward<T>(t));
}

template <typename T>
auto operator,(T &t, wrap<void, false>) {
    static_assert(!std::is_reference<T>::value, "TODO Add message");
    DEBUG("T& called");
    return wrap<T, true>(t);
}

template <typename F>
auto call(F &&f) -> decltype(std::forward<F>(f)()) {
    auto result = (std::forward<F>(f)(), wrap_void());
    DEBUG("Function called");
    return result.get();
}

}  // namespace detail

inline bool init() {
    int provided;
    MPI_Init_thread(nullptr, nullptr, MPI_THREAD_FUNNELED, &provided);
    return provided >= MPI_THREAD_FUNNELED;
}

inline bool finalize() { return MPI_Finalize() == MPI_SUCCESS; }

struct MsgInfo {
    size_t count;
    bool flag;
};

// TODO make partask_mpi.cpp
// buffers should have sizeof(MsgInfo) free bytes before the beginning!
inline void mpi_send_buffer(char *buffer, size_t count, int destination, int tag, bool flag) {
    ASSERT_MAIN_THREAD;
    TIME_TRACE_SCOPE("partask::mpi_send_buffer");
    DEBUG("mpi_send_buffer() called");
    MsgInfo info{count, flag};
    memcpy(buffer - sizeof(info), &info, sizeof(info));
    size_t all_count = count + sizeof(info);
    VERIFY(all_count <= std::numeric_limits<int>::max());
    int rc = MPI_Send(buffer - sizeof(info), static_cast<int>(all_count), MPI_BYTE, destination, tag, MPI_COMM_WORLD);
    VERIFY(rc == MPI_SUCCESS);
}

inline MsgInfo mpi_recv_buffer(char *buffer, size_t buffer_size, int source, int tag) {
    DEBUG("mpi_recv_buffer() called");
    TIME_TRACE_SCOPE("partask::mpi_recv_buffer");
    size_t all_count = buffer_size + sizeof(MsgInfo);
    VERIFY(all_count <= std::numeric_limits<int>::max());
    MPI_Status status;
    int rc = MPI_Recv(buffer - sizeof(MsgInfo), static_cast<int>(all_count), MPI_BYTE, source, tag, MPI_COMM_WORLD, &status);
    VERIFY(rc == MPI_SUCCESS);
    int actual_count;
    MPI_Get_count(&status, MPI_BYTE, &actual_count);
    size_t count = actual_count - sizeof(MsgInfo);
    VERIFY(count <= buffer_size);
    MsgInfo info;
    memcpy(&info, buffer - sizeof(info), sizeof(info));
    return info;
}

inline void mpi_send_buffer_bcast(char *buffer, size_t count, size_t buffer_size, int root, bool flag) {
    ASSERT_MAIN_THREAD;
    TIME_TRACE_SCOPE("partask::mpi_send_buffer_bcast");
    DEBUG("mpi_send_buffer_bcast() called");
    MsgInfo info{count, flag};
    VERIFY(info.count || info.flag);
    memcpy(buffer - sizeof(info), &info, sizeof(info));
    size_t all_count = buffer_size + sizeof(info);
    VERIFY(all_count <= std::numeric_limits<int>::max());
    int rc = MPI_Bcast(buffer - sizeof(info), static_cast<int>(all_count), MPI_BYTE, root, MPI_COMM_WORLD);
    VERIFY(rc == MPI_SUCCESS);
}

inline MsgInfo mpi_recv_buffer_bcast(char *buffer, size_t buffer_size, int root) {
    ASSERT_MAIN_THREAD;
    TIME_TRACE_SCOPE("partask::mpi_recv_buffer_bcast");
    DEBUG("mpi_recv_buffer_bcast() called");
    size_t all_count = buffer_size + sizeof(MsgInfo);
    int rc = MPI_Bcast(buffer - sizeof(MsgInfo), static_cast<int>(all_count), MPI_BYTE, root, MPI_COMM_WORLD);  // count should be the same!
    VERIFY(rc == MPI_SUCCESS);
    MsgInfo info;
    memcpy(&info, buffer - sizeof(MsgInfo), sizeof(info));
    VERIFY(info.count || info.flag);
    return info;
}

// The example was taken from http://www.voidcn.com/article/p-vjnlygmc-gy.html
// Probably we should send/recv async
template <bool BROADCAST>
class OutputMPIBuffer : public std::streambuf {
    static const size_t KiloByte = 1 << 10;
    static const size_t MegaByte = 1 << 20;
    static const size_t GigaByte = 1 << 30;

public:
    // the buffer size should be EXACTELY the same as in correspondent InputBuffer
    explicit OutputMPIBuffer(int destination, int tag = 0, size_t buffer_size = 100 * MegaByte)
        : buffer_size_{buffer_size}, buffer_(buffer_size + sizeof(MsgInfo)), destination_{destination}, tag_{tag} {
        VERIFY(buffer_.size() <= std::numeric_limits<int>::max());
        setp(data(), data() + buffer_size_);
    }

    ~OutputMPIBuffer() noexcept {
        flush(true);
        // INFO("Flush called: " << flush_count_);
    }
    OutputMPIBuffer(OutputMPIBuffer &&) = default;
    OutputMPIBuffer &operator=(OutputMPIBuffer &&) = default;

protected:
    int sync() override {
        flush(false);
        return 0;
    }

    int_type overflow(int_type ch) override {
        flush(false);
        if (ch != traits_type::eof()) {
            VERIFY(std::less_equal<char *>()(pptr(), epptr()));
            *pptr() = traits_type::to_char_type(ch);
            pbump(1);
        }
        return ch;
    }

private:
    char *data() {
        return buffer_.data() + sizeof(MsgInfo);
    }

    void resize_buffer(size_t buffer_size) {
        buffer_size_ = buffer_size;
        buffer_.resize(buffer_size_ + sizeof(MsgInfo));
    }

    void flush(bool last) {
        // ++flush_count_;
        size_t count = pptr() - pbase();
        if (!count && !last) {
            return;
        }
        if (BROADCAST) {
            mpi_send_buffer_bcast(pbase(), count, buffer_size_, destination_, last);
        } else {
            mpi_send_buffer(pbase(), count, destination_, tag_, last);
        }
        setp(data(), data() + buffer_size_);
    }

    OutputMPIBuffer(const OutputMPIBuffer &) = delete;
    OutputMPIBuffer &operator=(const OutputMPIBuffer &) = delete;

private:
    size_t buffer_size_;
    std::vector<char> buffer_;
    int destination_, tag_;
    // size_t flush_count_ = 0;
};

template <bool BROADCAST>
class InputMPIBuffer : public std::streambuf {
    static const size_t KiloByte = 1 << 10;
    static const size_t MegaByte = 1 << 20;
    static const size_t GigaByte = 1 << 30;

public:
    // the buffer size should be EXACTELY the same as in correspondent OutputBuffer
    explicit InputMPIBuffer(int source, int tag = MPI_ANY_TAG, size_t buffer_size = 100 * MegaByte, size_t put_back = 1 * KiloByte)
        : buffer_size_{buffer_size}, put_back_{std::max(put_back, sizeof(MsgInfo))}, buffer_(put_back_ + buffer_size),
          source_{source}, tag_{tag}, last_{false} {
        VERIFY(buffer_.size() <= std::numeric_limits<int>::max());
        setg(buffer_.data(), data() + buffer_size_, data() + buffer_size_);
    }

    ~InputMPIBuffer() noexcept {
        while (pull()) {}
    }

    InputMPIBuffer(InputMPIBuffer &&) = default;
    InputMPIBuffer &operator=(InputMPIBuffer &&) = default;

private:
    char *data() {
        return buffer_.data() + put_back_;
    }

    void resize_buffer(size_t buffer_size) {
        buffer_size_ = buffer_size;
        buffer_.resize(buffer_size_ + put_back_);
    }

    MsgInfo retrive() {
        MsgInfo info;
        if (BROADCAST) {
            info = mpi_recv_buffer_bcast(data(), buffer_size_, source_);
        } else {
            info = mpi_recv_buffer(data(), buffer_size_, source_, tag_);
        }
        return info;
    }

    size_t pull() {
        if (last_) {
            return 0;
        }

        MsgInfo info;
        do {
            info = retrive();
            VERIFY(info.flag || info.count);
        } while (!info.flag && info.count == 0);

        last_ = info.flag;
        return info.count;
    }


    int_type underflow() override {
        if (gptr() < egptr()) {  // buffer not exhausted
            return traits_type::to_int_type(*gptr());
        }

        size_t n = pull();
        if (n == 0) return traits_type::eof();

        // Set buffer pointers
        setg(buffer_.data(), data(), data() + n);

        return traits_type::to_int_type(*gptr());
    }

    InputMPIBuffer(const InputMPIBuffer &) = delete;
    InputMPIBuffer &operator=(const InputMPIBuffer &) = delete;
private:
    size_t buffer_size_;
    size_t put_back_;
    std::vector<char> buffer_;
    int source_, tag_;
    bool last_;
};

inline void mpi_send_buffer_async(char *buffer, size_t count, int destination, int tag, bool flag, MPI_Request &req) {
    ASSERT_MAIN_THREAD;
    TIME_TRACE_SCOPE("partask::mpi_send_buffer_async");
    DEBUG("mpi_send_buffer() called");
    MsgInfo info{count, flag};
    memcpy(buffer - sizeof(info), &info, sizeof(info));
    size_t all_count = count + sizeof(info);
    VERIFY(all_count <= std::numeric_limits<int>::max());
    int rc = MPI_Isend(buffer - sizeof(info), static_cast<int>(all_count), MPI_BYTE, destination, tag, MPI_COMM_WORLD, &req);
    VERIFY(rc == MPI_SUCCESS);
}

inline void mpi_recv_buffer_async(char *buffer, size_t buffer_size, int source, int tag, MPI_Request &req) {
    DEBUG("mpi_recv_buffer() called");
    TIME_TRACE_SCOPE("partask::mpi_recv_buffer_async");
    size_t all_count = buffer_size + sizeof(MsgInfo);
    VERIFY(all_count <= std::numeric_limits<int>::max());
    int rc = MPI_Irecv(buffer - sizeof(MsgInfo), static_cast<int>(all_count), MPI_BYTE, source, tag, MPI_COMM_WORLD, &req);
    VERIFY(rc == MPI_SUCCESS);
}

inline MsgInfo mpi_recv_buffer_wait(char *buffer, MPI_Request &req) {
    TIME_TRACE_SCOPE("partask::mpi_recv_buffer_wait");
    MPI_Status status;
    MPI_Wait(&req, &status);
    int actual_count;
    MPI_Get_count(&status, MPI_BYTE, &actual_count);
    size_t count = actual_count - sizeof(MsgInfo);
    MsgInfo info;
    memcpy(&info, buffer - sizeof(info), sizeof(info));
    VERIFY(count == info.count);
    return info;
}

inline void mpi_send_buffer_bcast_async(char *buffer, size_t count, size_t buffer_size, int root, bool flag, MPI_Request &req) {
    ASSERT_MAIN_THREAD;
    TIME_TRACE_SCOPE("partask::mpi_send_buffer_bcast_async");
    DEBUG("mpi_send_buffer_bcast_async() called. count = " << count << " flag " << flag);
    MsgInfo info{count, flag};
    VERIFY(info.count || info.flag);
    memcpy(buffer - sizeof(info), &info, sizeof(info));
    size_t all_count = buffer_size + sizeof(info);
    VERIFY(all_count <= std::numeric_limits<int>::max());
    int rc = MPI_Ibcast(buffer - sizeof(info), static_cast<int>(all_count), MPI_BYTE, root, MPI_COMM_WORLD, &req);
    VERIFY(rc == MPI_SUCCESS);
}

inline void mpi_recv_buffer_bcast_async(char *buffer, size_t buffer_size, int root, MPI_Request &req) {
    ASSERT_MAIN_THREAD;
    TIME_TRACE_SCOPE("partask::mpi_recv_buffer_bcast_async");
    DEBUG("mpi_recv_buffer_bcast() called");
    size_t all_count = buffer_size + sizeof(MsgInfo);
    int rc = MPI_Ibcast(buffer - sizeof(MsgInfo), static_cast<int>(all_count), MPI_BYTE, root, MPI_COMM_WORLD, &req);  // count should be the same!
    VERIFY(rc == MPI_SUCCESS);
}

inline MsgInfo mpi_recv_buffer_bcast_wait(char *buffer, MPI_Request &req) {
    TIME_TRACE_SCOPE("partask::mpi_recv_buffer_bcast_wait");
    MPI_Wait(&req, MPI_STATUS_IGNORE);
    MsgInfo info;
    memcpy(&info, buffer - sizeof(MsgInfo), sizeof(info));
    DEBUG("mpi_recv_buffer_bcast_wait() called. count = " << info.count << " flag " << info.flag);
    VERIFY(info.count || info.flag);
    return info;
}

template <bool BROADCAST>
class OutputMPIBufferAsync : public std::streambuf {
    static const size_t KiloByte = 1 << 10;
    static const size_t MegaByte = 1 << 20;
    static const size_t GigaByte = 1 << 30;

public:
    // the buffer size should be EXACTELY the same as in correspondent InputBuffer
    explicit OutputMPIBufferAsync(size_t destination, int tag = 0, size_t buffer_size = 100 * MegaByte)
        : buffer_size_{buffer_size}, buffer1_(buffer_size + sizeof(MsgInfo)), buffer2_(buffer_size + sizeof(MsgInfo)), destination_{int(destination)}, tag_{tag} {
        pbuffer_ = &buffer1_;
        setp(data(), data() + buffer_size_);
    }

    ~OutputMPIBufferAsync() noexcept {
        flush(true);
        if (wait_) {
            MPI_Wait(&req_, MPI_STATUS_IGNORE);
            wait_ = false;
        }
        // INFO("Flush called: " << flush_count_);
    }
    OutputMPIBufferAsync(OutputMPIBufferAsync &&) = default;
    OutputMPIBufferAsync &operator=(OutputMPIBufferAsync &&) = default;

protected:
    int sync() override {
        flush(false);
        return 0;
    }

    int_type overflow(int_type ch) override {
        flush(false);
        if (ch != traits_type::eof()) {
            VERIFY(std::less_equal<char *>()(pptr(), epptr()));
            *pptr() = traits_type::to_char_type(ch);
            pbump(1);
        }
        return ch;
    }

private:
    char *data() {
        return pbuffer_->data() + sizeof(MsgInfo);
    }

    void flush(bool last) {
        // ++flush_count_;
        size_t count = pptr() - pbase();
        if (!count && !last) {
            return;
        }
        if (wait_) {
            MPI_Wait(&req_, MPI_STATUS_IGNORE);
            wait_ = false;
        }
        if (BROADCAST) {
            mpi_send_buffer_bcast_async(pbase(), count, buffer_size_, destination_, last, req_);
        } else {
            mpi_send_buffer_async(pbase(), count, destination_, tag_, last, req_);
        }
        wait_ = true;
        pbuffer_ = pbuffer_ == &buffer1_ ? &buffer2_ : &buffer1_;
        setp(data(), data() + buffer_size_);
    }

    OutputMPIBufferAsync(const OutputMPIBufferAsync &) = delete;
    OutputMPIBufferAsync &operator=(const OutputMPIBufferAsync &) = delete;

private:
    size_t buffer_size_;
    std::vector<char> buffer1_, buffer2_;
    int destination_, tag_;
    // size_t flush_count_ = 0;
    std::vector<char> *pbuffer_;
    bool wait_ = false;
    MPI_Request req_;
};

template <bool BROADCAST>
class InputMPIBufferAsync : public std::streambuf {
    static const size_t KiloByte = 1 << 10;
    static const size_t MegaByte = 1 << 20;
    static const size_t GigaByte = 1 << 30;

public:
    // the buffer size should be EXACTELY the same as in correspondent OutputBuffer
    explicit InputMPIBufferAsync(size_t source, int tag = MPI_ANY_TAG, size_t buffer_size = 100 * MegaByte, size_t put_back = 1 * KiloByte)
        : buffer_size_{buffer_size}, put_back_{std::max(put_back, sizeof(MsgInfo))}, buffer1_(put_back_ + buffer_size), buffer2_(put_back_ + buffer_size),
          source_{int(source)}, tag_{tag} {
        setg(pbuffer_->data(), data() + buffer_size_, data() + buffer_size_);
    }

    ~InputMPIBufferAsync() noexcept {
        while (pull()) {}
    }

    InputMPIBufferAsync(InputMPIBufferAsync &&) = default;
    InputMPIBufferAsync &operator=(InputMPIBufferAsync &&) = default;

protected:
    int_type underflow() override {
        if (gptr() < egptr()) {  // buffer not exhausted
            return traits_type::to_int_type(*gptr());
        }

        size_t n = pull();
        if (n == 0) return traits_type::eof();

        // Set buffer pointers
        setg(pbuffer_->data(), data(), data() + n);

        return traits_type::to_int_type(*gptr());
    }

private:
    char *data() {
        return pbuffer_->data() + put_back_;
    }

    MsgInfo wait() {
        return BROADCAST ? mpi_recv_buffer_bcast_wait(data(), req_) : mpi_recv_buffer_wait(data(), req_);
    }

    void recv() {
        if (BROADCAST) {
            mpi_recv_buffer_bcast_async(data(), buffer_size_, source_, req_);
        } else {
            mpi_recv_buffer_async(data(), buffer_size_, source_, tag_, req_);
        }
    }

    size_t pull() {
        if (last_) {
            return 0;
        }

        pbuffer_ = pbuffer_ == &buffer1_ ? &buffer2_ : &buffer1_;
        if (!wait_) {
            recv();
            wait_ = true;
        }

        VERIFY(wait_);
        MsgInfo info = wait();
        last_ = info.flag;
        wait_ = false;

        if (!last_) {
            pbuffer_ = pbuffer_ == &buffer1_ ? &buffer2_ : &buffer1_;
            recv();
            pbuffer_ = pbuffer_ == &buffer1_ ? &buffer2_ : &buffer1_;
            wait_ = true;
        }

        return info.count;
    }

    InputMPIBufferAsync(const InputMPIBufferAsync &) = delete;
    InputMPIBufferAsync &operator=(const InputMPIBufferAsync &) = delete;

private:
    size_t buffer_size_;
    size_t put_back_;
    std::vector<char> buffer1_, buffer2_;
    int source_, tag_;
    bool last_ = false;
    bool wait_ = false;
    std::vector<char> *pbuffer_ = &buffer1_;
    MPI_Request req_;
};

// TODO add put_back support
// TODO Add mutex
class ChunkedStringBuffer : public std::streambuf {
    static const size_t KiloByte = 1 << 10;
    static const size_t MegaByte = 1 << 20;
    static const size_t GigaByte = 1 << 30;

public:
    explicit ChunkedStringBuffer(size_t buffer_size = 100 * MegaByte)
        : buffer_size_{buffer_size}, g_buffer_id_{size_t(-1)}, buffer_(buffer_size_) {
        setp(buffer_.data(), buffer_.data() + buffer_size_);
        setg(buffer_.data(), buffer_.data() + buffer_size_, buffer_.data() + buffer_size_);
    }

    int sync() override {
        flush();
        return 0;
    }

    size_t size() const {
        size_t result = pptr() - pbase();
        for (const auto &buffer : buffers_) {
            result += buffer.size();
        }
        return result;
    }

    int_type underflow() override {
        if (gptr() < egptr()) {  // buffer not exhausted
            return traits_type::to_int_type(*gptr());
        }

        /* DEBUG("underflow(): g_buffer_id_ = " << g_buffer_id_); */
        if (g_buffer_id_ != size_t(-1)) {
            // clear current buffer
            buffers_[g_buffer_id_].clear();
            buffers_[g_buffer_id_].shrink_to_fit();
        }

        ++g_buffer_id_;

        if (g_buffer_id_ == buffers_.size()) {
            flush();
            VERIFY(g_buffer_id_ < buffers_.size());
            if (buffers_[g_buffer_id_].empty()) {
                return traits_type::eof();
            }
        }

        // Set buffer pointers
        auto &b = buffers_[g_buffer_id_];
        setg(b.data(), b.data(), b.data() + b.size());

        return traits_type::to_int_type(*gptr());
    }

    int_type overflow(int_type ch) override {
        flush();

        if (ch != traits_type::eof()) {
            VERIFY(std::less_equal<char *>()(pptr(), epptr()));
            *pptr() = traits_type::to_char_type(ch);
            pbump(1);
        }
        return ch;
    }

    void broadcast(int root = 0) {
        if (world_rank() == root) {
            flush();
        }
        size_t buffers_size = buffers_.size();
        ::partask::broadcast(buffers_size, root);
        buffers_.resize(buffers_size);

        std::vector<size_t> sizes(buffers_size + 1 + 1 + 1 + 6);
        sizes[0] = buffer_size_;  // Actually already should be equal
        sizes[1] = g_buffer_id_;
        sizes[2] = buffer_.size();

        sizes[3] = eback() - buffer_.data();
        sizes[4] = gptr() - buffer_.data();
        sizes[5] = egptr() - buffer_.data();
        sizes[6] = pbase() - buffer_.data();
        sizes[7] = pptr() - buffer_.data();  // Actually, could not be synced
        sizes[8] = epptr() - buffer_.data();
        for (size_t i = 0; i < buffers_.size(); ++i) {
            sizes[i + 9] = buffers_[i].size();
        }
        ::partask::membroadcast(sizes.data(), sizeof(sizes[0]) * sizes.size(), root);
        buffer_size_ = sizes[0];
        g_buffer_id_ = sizes[1];
        buffer_.resize(sizes[2]);
        ::partask::membroadcast(buffer_.data(), sizeof(buffer_[0]) * buffer_.size(), root);
        if (world_rank() != root) {
            setg(sizes[3] + buffer_.data(), sizes[4] + buffer_.data(), sizes[5] + buffer_.data());
            setp(sizes[6] + buffer_.data(), sizes[8] + buffer_.data());
        }

        for (size_t i = 0; i < buffers_.size(); ++i) {
            buffers_[i].resize(sizes[i + 9]);
            ::partask::membroadcast(buffers_[i].data(), sizeof(buffers_[i][0]) * buffers_[i].size(), root);
        }
    }

private:
    void flush() {
        size_t count = pptr() - pbase();
        buffer_.resize(count);
        buffers_.emplace_back(buffer_size_);
        std::swap(buffer_, buffers_.back());
        setp(buffer_.data(), buffer_.data() + buffer_size_);
    }

    size_t buffer_size_;
    size_t g_buffer_id_;
    std::vector<std::vector<char>> buffers_;
    std::vector<char> buffer_;
};

template <typename BaseStream, typename Buffer>
class MPIStream : public BaseStream {
    using This = MPIStream<BaseStream, Buffer>;
public:
    template <typename... Args>
    MPIStream(Args... args) : BaseStream(nullptr), buf_(args...) {  // All args are integral
        this->init(&buf_);
    }

    MPIStream(This&&) = default;
    MPIStream& operator=(This&&) = default;
    MPIStream(const This&) = delete;
    MPIStream& operator=(const This&) = delete;

private:
    Buffer buf_;
};

using InputMPIStream = MPIStream<std::istream, InputMPIBufferAsync<false>>;
using OutputMPIStream = MPIStream<std::ostream, OutputMPIBufferAsync<false>>;
class InputMPIStreamBcast : public MPIStream<std::istream, InputMPIBufferAsync<true>> {
public:
    using MPIStream::MPIStream;
    InputMPIStreamBcast() : MPIStream(0) {}
};

class OutputMPIStreamBcast : public MPIStream<std::ostream, OutputMPIBufferAsync<true>> {
public:
    using MPIStream::MPIStream;
    OutputMPIStreamBcast() : MPIStream(0) {}
};

class ChunkedStringStream : public MPIStream<std::iostream, ChunkedStringBuffer> {
public:
    size_t size() const { return dynamic_cast<ChunkedStringBuffer *>(this->rdbuf())->size(); }

    void broadcast(int root = 0) { return dynamic_cast<ChunkedStringBuffer *>(this->rdbuf())->broadcast(root); }
};

template <typename T, typename Serialize, typename Deserialize>
void broadcast(T &data, Serialize &&serialize, Deserialize &&deserialize, int root = 0) {
    TIME_TRACE_SCOPE("partask::broadcast", utils::type_name<T>());

    ASSERT_MAIN_THREAD;
    DEBUG("Broadcasting of type " << typeid(T).name());

    static size_t call_count = 0;
    DEBUG("membroadcast() called " << call_count << " times");
    ++call_count;

    if (world_rank() == root) {
        OutputMPIStreamBcast os(root);
        DEBUG("Broadcast serialization...");
        std::forward<Serialize>(serialize)(os, data);
        DEBUG("Broadcast serialization complete");
    } else {
        InputMPIStreamBcast is(root);
        DEBUG("Broadcast deserialization...");
        std::forward<Deserialize>(deserialize)(is, data);
        DEBUG("Broadcast deserialization complete");
    }
}

template <typename T1, typename Serialize1, typename Deserialize1, typename T2, typename Serialize2, typename Deserialize2>
void broadcast2(T1 &data1, Serialize1 &&serialize1, Deserialize1 &&deserialize1,
                T2 &data2, Serialize2 &&serialize2, Deserialize2 &&deserialize2,
                int root = 0) {
    TIME_TRACE_SCOPE("partask::broadcast2", utils::type_name<T1>() + "+" + utils::type_name<T2>());

    ASSERT_MAIN_THREAD;
    DEBUG("Broadcasting of types " << typeid(T1).name() << " " << typeid(T2).name());

    static size_t call_count = 0;
    DEBUG("membroadcast() called " << call_count << " times");
    ++call_count;

    if (world_rank() == root) {
        OutputMPIStreamBcast os(root);
        DEBUG("Broadcast serialization...");
        std::forward<Serialize1>(serialize1)(os, data1);
        std::forward<Serialize2>(serialize2)(os, data2);
        DEBUG("Broadcast serialization complete");
    } else {
        InputMPIStreamBcast is(root);
        DEBUG("Broadcast deserialization...");
        std::forward<Deserialize1>(deserialize1)(is, data1);
        std::forward<Deserialize2>(deserialize2)(is, data2);
        DEBUG("Broadcast deserialization complete");
    }
}

template <typename T, typename Serialize, typename Deserialize>
void broadcast_full_dump(T &data, Serialize &&serialize, Deserialize &&deserialize, int root = 0) {
    ASSERT_MAIN_THREAD;
    TIME_TRACE_SCOPE("partask::broadcast_full_dump", utils::type_name<T>());

    DEBUG("Broadcasting of type " << typeid(T).name());

    static size_t call_count = 0;
    DEBUG("broadcast_full_dump() called " << call_count << " times");
    ++call_count;

    ChunkedStringStream css;
    if (world_rank() == root) {
        DEBUG("Broadcast serialization...");
        std::forward<Serialize>(serialize)(css, data);
        DEBUG("Broadcast serialization complete");
    }
    css.broadcast(root);

    if (world_rank() != root) {
        DEBUG("Broadcast deserialization...");
        std::forward<Deserialize>(deserialize)(css, data);
        DEBUG("Broadcast deserialization complete");
    }
}

template <typename T>
auto broadcast(T &data, int root = 0) -> decltype(std::declval<std::enable_if_t<!std::is_pod<T>::value, T>>(),
                                                  io::binary::BinWrite(declref<std::ostream>(), data),
                                                  io::binary::BinRead(declref<std::istream>(), data),
                                                  void()) {
    broadcast(data, [](std::ostream &os, const T &data) { io::binary::BinWrite(os, data); },
              [](std::istream &is, T &data) { io::binary::BinRead(is, data); },
              root);
}

template <typename T, typename Serialize>
auto send(const T &data, Serialize &&serialize, int destination, int tag = 0) -> decltype(std::forward<Serialize>(serialize)(declref<std::ostream>(), data), void()) {
    ASSERT_MAIN_THREAD;
    TIME_TRACE_SCOPE("partask::send", utils::type_name<T>());

    OutputMPIStream os(destination, tag);
    DEBUG("Serialization...");
    std::forward<Serialize>(serialize)(os, data);
}

template <typename T, typename Deserialize>
auto recv(T &data, Deserialize &&deserialize, int source, int tag = MPI_ANY_TAG) -> decltype(std::forward<Deserialize>(deserialize)(declref<std::istream>(), data), void()) {
    ASSERT_MAIN_THREAD;
    TIME_TRACE_SCOPE("partask::recv", utils::type_name<T>());

    InputMPIStream is(source, tag);
    DEBUG("Serialization...");
    std::forward<Deserialize>(deserialize)(is, data);
}

inline std::vector<int> collect_num_threads(int root = 0) {
    ASSERT_MAIN_THREAD;
    std::vector<int> all_num_threads;

    int num_threads = omp_get_max_threads();
    if (world_rank() == root) {
        all_num_threads.resize(world_size());
        MPI_Gather(&num_threads, 1, MPI_INT, all_num_threads.data(), 1, MPI_INT, root, MPI_COMM_WORLD);
    } else {
        MPI_Gather(&num_threads, 1, MPI_INT, nullptr, 1, MPI_INT, root, MPI_COMM_WORLD);
    }
    broadcast(all_num_threads, root);
    return all_num_threads;
}

inline int overall_num_threads(int root = 0) {
    ASSERT_MAIN_THREAD;
    auto threads = collect_num_threads(root);
    return std::accumulate(threads.cbegin(), threads.cend(), int(0));
}

inline void all_set_num_threads(const std::vector<int> &all_num_threads, int root = 0) {
    ASSERT_MAIN_THREAD;
    int num_threads;
    MPI_Scatter(const_cast<int *>(all_num_threads.data()), 1, MPI_INT, &num_threads, 1, MPI_INT, root, MPI_COMM_WORLD);
    omp_set_num_threads(num_threads);
}

inline void all_set_num_threads(int num_threads, int = 0) { omp_set_num_threads(num_threads); }

#define CREATE_HAS_METHOD_CHECKER(METHOD) \
struct has_##METHOD##_method { \
    template <typename T, typename Tuple, std::size_t... I> \
    constexpr static auto detail_test(T &&obj, Tuple &&tuple, std::index_sequence<I...>) -> \
        decltype(std::forward<T>(obj).METHOD(std::get<I>(std::forward<Tuple>(tuple))...), \
                    std::true_type()); \
\
    template <typename T, typename Tuple> \
    constexpr static auto test(T &&obj, Tuple &&tuple) -> decltype(detail_test(std::forward<T>(obj), \
                                                                                std::forward<Tuple>(tuple),\
                                                                                std::make_index_sequence<std::tuple_size_v<std::remove_reference_t<Tuple>>>{})); \
    constexpr static std::false_type test(...); \
};

class TaskRegistry {
    static const int MAP_TAG = 13;
    static const int MERGE_TAG = 14;
    static const size_t MULT = 1;
    static const size_t STOP_LISTENING_TASK = -1;

    class AbstractTask {
    public:
        virtual void process(std::istream &is, std::ostream &) {
            // Skip all characters in stream
            DEBUG("Unimplemented (trivial) process method");
            is.ignore(std::numeric_limits<std::streamsize>::max());
        }

        virtual void sync(void) {
            DEBUG("Unimplemented (trivial) sync method");
            // Do nothing
        }

        virtual ~AbstractTask() noexcept = default;
    };

    CREATE_HAS_METHOD_CHECKER(sync);
    CREATE_HAS_METHOD_CHECKER(merge);
    CREATE_HAS_METHOD_CHECKER(process);
    CREATE_HAS_METHOD_CHECKER(make_splitter);

    class AbstractTaskFactory {
    public:
        virtual AbstractTask *create(std::istream&) const = 0;
        virtual ~AbstractTaskFactory() noexcept = default;
    };

    template <typename Task, typename Locals>
    class TaskFactory : public AbstractTaskFactory {
        class ConcreteTask : public AbstractTask {
            constexpr static const bool has_make_splitter = decltype(has_make_splitter_method::test(std::declval<Task>(),
                                                                                                    std::tuple_cat(std::tuple<size_t>(),
                                                                                                                   std::declval<Locals>())))::value;
            constexpr static const bool has_process = decltype(has_process_method::test(std::declval<Task>(),
                                                                                        std::tuple_cat(std::tuple<std::istream&, std::ostream&>(declref<std::istream>(),
                                                                                                                                                declref<std::ostream>()),
                                                                                                       std::declval<Locals>())))::value;
            constexpr static const bool has_merge = decltype(has_merge_method::test(std::declval<Task>(),
                                                                                    std::tuple_cat(std::tuple<std::vector<std::istream*>>(),
                                                                                                   std::declval<Locals>())))::value;
            constexpr static const bool has_sync = decltype(has_sync_method::test(std::declval<Task>(), std::declval<Locals>()))::value;
        public:
            ConcreteTask(Task &&task, const Locals &locals) : task_{std::move(task)}, locals_{locals} {}

            template <typename SIZE = size_t>
            auto make_splitter(std::enable_if_t<has_make_splitter, SIZE> size) {
                auto make_splitter_args = std::tuple_cat(std::make_tuple(size), locals_);
                auto splitter =
                        std::apply(
                            [&](auto &&... ts) {
                                return task_.make_splitter(std::forward<decltype(ts)>(ts)...);
                            },
                            make_splitter_args);
                return splitter;
            }

            template <typename SIZE = size_t>
            auto make_splitter(std::enable_if_t<!has_make_splitter, SIZE>) {
                WARN("Unimplemented (trivial) make_splitter method");
                auto trivial_splitter = [](std::ostream&, int) {
                    return false;
                };
                return trivial_splitter;
            }

            void process(std::istream &is, std::ostream &os) override {
                TIME_TRACE_SCOPE("partask::task:process");
                TIME_TRACE_SCOPE(utils::type_name<Task>());
                process_impl(is, os);
            }

            void sync(void) override {
                TIME_TRACE_SCOPE("partask::task:sync");
                TIME_TRACE_SCOPE(utils::type_name<Task>());
                sync_impl();
            }

            template <typename T = const std::vector<std::istream*>>
            decltype(auto) merge(std::enable_if_t<has_merge, T> &piss) {
                TIME_TRACE_SCOPE("partask::task:merge");
                TIME_TRACE_SCOPE(utils::type_name<Task>());
                auto merge_args = std::tuple_cat(std::make_tuple(piss), locals_);
                auto merge_call = [&](auto &&... ts) { return task_.merge(std::forward<decltype(ts)>(ts)...); };
                return std::apply(merge_call, merge_args);
            }

            template <typename T = const std::vector<std::istream *>>
            void merge(std::enable_if_t<!has_merge, T> &piss) {
                DEBUG("Unimplemented (trivial) merge method");
                // Do nothing except reading streams
                for (auto &pis : piss) {
                    auto &is = *pis;
                    // Skip all characters in stream
                    is.ignore(std::numeric_limits<std::streamsize>::max());
                }
            }

        private:
            Task task_;
            Locals locals_;

            template <typename ISTREAM = std::istream>
            void process_impl(std::enable_if_t<has_process, ISTREAM> &is, std::ostream &os) {
                auto args = std::tuple_cat(std::tuple<std::istream &, std::ostream &>(is, os), locals_);
                auto process_call = [this](auto &&... ts) { return task_.process(std::forward<decltype(ts)>(ts)...); };
                std::apply(process_call, args);
            }

            template <typename ISTREAM = std::istream>
            void process_impl(std::enable_if_t<!has_process, ISTREAM> &is, std::ostream &os) {
                AbstractTask::process(is, os);
            }

            template <typename VOID = void>
            std::enable_if_t<has_sync, VOID> sync_impl() {
                TIME_TRACE_SCOPE(utils::type_name<Task>());
                auto sync_call = [this](auto &&... ts) { return task_.sync(std::forward<decltype(ts)>(ts)...); };
                std::apply(sync_call, locals_);
            }

            template <typename VOID = void>
            std::enable_if_t<!has_sync, VOID> sync_impl() {
                AbstractTask::sync();
            }

        };

    public:
        TaskFactory(Locals &&locals) : locals_{std::move(locals)} {}

        ConcreteTask *acquire(Task &&task) const {
            return new ConcreteTask(std::move(task), locals_);
        }

        ConcreteTask *create(std::istream& is) const override {
            Task task(is);
            return acquire(std::move(task));
        }

    private:
        Locals locals_;
    };


public:
    TaskRegistry() {
        ASSERT_MAIN_THREAD;
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank_);
        MPI_Comm_size(MPI_COMM_WORLD, &world_size_);
        all_num_threads_ = collect_num_threads();
    }

    template <typename Task, typename Locals>
    class Job {
    public:
        Job(TaskRegistry &task_registry, size_t job_id)
                : task_registry_{task_registry}, job_id_{job_id} {}

        template <typename... Args>
        decltype(auto) operator()(Args &&... args) const {
            TIME_TRACE_SCOPE("partask::job");
            VERIFY(task_registry_.world_rank_ == 0);
            Task task(std::forward<Args>(args)...);

            task_registry_.job_broadcast_(job_id_);
            {
                OutputMPIStreamBcast obs(0);
                task.serialize(obs);
            } // close obs stream and send the rest of the data

            int world_size;
            MPI_Comm_size(MPI_COMM_WORLD, &world_size);

            auto pfactory = dynamic_cast<TaskFactory<Task, Locals>*>(task_registry_.factories_[job_id_].get());
            auto ptask = as_unique(pfactory->acquire(std::move(task)));

            auto plocal_map = std::make_shared<ChunkedStringStream>();
            auto plocal_merge = std::make_shared<ChunkedStringStream>();
            {
                std::vector<std::shared_ptr<std::ostream>> oss;
                oss.push_back(plocal_map);
                for (int rank = 1; rank < world_size; ++rank) {
                    oss.emplace_back(new OutputMPIStream(rank, MAP_TAG));
                }

                const auto &all_num_threads_ = task_registry_.all_num_threads_;
                size_t sum_num_threads = std::accumulate(all_num_threads_.cbegin(), all_num_threads_.cend(), 0);
                VERIFY(sum_num_threads > 0);
                DEBUG("All threads: " << sum_num_threads << " Multiplicator: " << MULT);
                auto splitter = ptask->make_splitter(sum_num_threads * MULT);

                auto mult_splitter = [&splitter](auto &os, int rank, size_t count) {
                    for (size_t i = 0; i < count; ++i) {
                        bool result = splitter(os, rank);
                        if (!result) return false;
                    }
                    return true;
                };

                for (int rank = 0; mult_splitter(*oss[rank], rank, all_num_threads_[rank]);
                     rank = (rank + 1) % world_size) {
                }
            }  // close streams here and send split data

            DEBUG("Process started");
            ptask->process(*plocal_map, *plocal_merge);
            DEBUG("Process done");

            std::vector<std::shared_ptr<std::istream>> iss;
            iss.push_back(std::move(plocal_merge));
            for (int rank = 1; rank < world_size; ++rank) {
                DEBUG("Getting data from node " << rank);
                iss.emplace_back(new InputMPIStream(rank, MERGE_TAG));
                DEBUG("Got data from node " << rank);
            }

            std::vector<std::istream *> piss;
            for (int rank = 0; rank < world_size; ++rank) {
                piss.push_back(iss[rank].get());
            }

            DEBUG("Merge started");
            auto wrap = (ptask->merge(piss), detail::wrap_void());
            DEBUG("Merge done, closing streams");
            for (auto is : iss) {
                is.reset();
            }
            DEBUG("Streams closed");
            DEBUG("sync calling...");
            ptask->sync();
            DEBUG("sync called, returning");
            return wrap.get();
        }

    private:
        TaskRegistry &task_registry_;
        size_t job_id_;
    };

    template <typename Task, typename... LocalArgs>
    auto add(LocalArgs &&... local_args) {
        barrier();

        auto locals = std::make_tuple(std::forward<LocalArgs>(local_args)...);  // std::make_tuple unwraps ref & cref
        size_t job_id = factories_.size();
        factories_.emplace_back(new TaskFactory<Task, decltype(locals)>(std::move(locals)));
        return Job<Task, decltype(locals)>(*this, job_id);
    }

    void stop_listening() {
        DEBUG("Stop listening");
        if (worker()) return;
        if (listening_) {
            job_broadcast_(STOP_LISTENING_TASK);
            listening_ = false;
        }
    }

    void stop() {
        WARN("stop() is depricated and will be removed soon. Use stop_listening() instead");
        return stop_listening();
    }

    void listen() {
        if (master()) {
            VERIFY(!listening_);
            listening_ = true;
            return;
        }
        DEBUG("Listening started");
        while (listen_one_()) {
        };
    }

    ~TaskRegistry() {
        stop_listening();
    }

    bool master() const { return world_rank() == 0; }
    bool worker() const { return !master(); }
    int world_size() const { return world_size_; }
    int world_rank() const { return world_rank_; }
    bool listening() const { return listening_; }

private:
    int world_rank_;
    int world_size_;
    std::vector<int> all_num_threads_;

    std::vector<std::unique_ptr<AbstractTaskFactory>> factories_;

    bool listening_ = false;

    void job_broadcast_(size_t job_id) {
        VERIFY(world_rank_ == 0);
        DEBUG("Job sending... " << job_id);
        broadcast(job_id);
        DEBUG("Job sent");
    }

    bool listen_one_() {
        VERIFY(worker());
        size_t job_id;
        DEBUG("Awaiting for job...");
        broadcast(job_id);  // Get job id
        DEBUG("Job got, id: " << job_id);
        if (job_id == STOP_LISTENING_TASK) {
            return false;
        }

        const auto &pfactory = factories_[job_id];
        DEBUG("Task object initializer obtained");
        auto pibs = std::make_unique<InputMPIStreamBcast>(0);
        auto ptask = as_unique(pfactory->create(*pibs));
        pibs.reset();  // Pull the rest of data and close the pipe
        DEBUG("Task object created");
        {
            InputMPIStream is(0, MAP_TAG);
            DEBUG("Input stream constructed");
            OutputMPIStream os(0, MERGE_TAG);
            DEBUG("Remote process output stream constructed");
            ptask->process(is, os);
            DEBUG("process done, closing the scope...");
        }  // close and destroy streams
        ptask->sync();
        DEBUG("sync() done");
        return true;
    }

};

#undef CREATE_HAS_METHOD_CHECKER

inline auto make_trivial_generator() {
    auto generator = [](std::ostream&, size_t) {
        return false;
    };
    return generator;
}

inline auto make_seq_generator(size_t size) {
    auto generator = [size, i = size_t(0)](std::ostream &os, size_t) mutable -> bool {
        if (i == size) {
            os.put('\0');
            return false;
        }
        os.put('\1');
        io::binary::BinWrite(os, i);
        ++i;
        return true;
    };

    return generator;
}

inline auto make_seq_plus_n_generator(size_t size) {
    auto generator = [size, i = size_t(0)](std::ostream &os, size_t) mutable -> bool {
        if (i == size) {
            return false;
        }
        io::binary::BinWrite(os, i, size);
        ++i;
        return true;
    };

    return generator;
}

inline std::vector<size_t> get_seq_plus_n(std::istream &is, size_t &size) {
    std::vector<size_t> chunks;
    size = 0;
    while (is.peek() != EOF) {
        size_t i;
        io::binary::BinRead(is, i, size);
        chunks.push_back(i);
    }

    return chunks;
}

template <typename Container>
auto make_seq_along_generator(const Container &c) {
    return make_seq_generator(c.size());
}

inline std::vector<size_t> get_seq(std::istream &is) {
    std::vector<size_t> chunks;
    // while (is.peek() != EOF) {
    while (is.get() && is) {
        size_t i;
        io::binary::BinRead(is, i);
        chunks.push_back(i);
    }

    return chunks;
}

template <typename T>
auto make_vector_splitter(size_t n, const std::vector<T>& data) {
    size_t N = data.size();

    auto splitter = [&data, n, N, chunk = size_t(0),
                     idx = size_t(0)](std::ostream& os, size_t /*node*/) mutable -> bool {
        if (chunk == n) {
            os.put('\0');
            return false;
        }

        size_t size = N / n + (chunk < (N % n));

        os.put('\1');
        io::binary::BinWrite(os, size);
        for (size_t i = 0; i < size; ++i, ++idx) {
            io::binary::BinWrite(os, data[idx]);
        }

        ++chunk;
        return true;
    };

    return splitter;
}

template <typename T>
auto all_equal(const T &v) -> decltype(broadcast(*new T(v)), T(v) == T(v), bool()) {
    TIME_TRACE_SCOPE("partask::all_equal");
    T cv(v);
    broadcast(cv);
    return v == cv;
}

template <typename F>
auto critical_ordered(F &&f) -> decltype(std::forward<F>(f)()) {
    TIME_TRACE_SCOPE("partask::critically_ordered");

    using wrap_type = decltype(std::forward<F>(f)(), detail::wrap_void());
    std::unique_ptr<wrap_type> pwrap;

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    for (int rank = 0; rank < world_size; ++rank) {
        barrier();

        if (world_rank == rank) {
            auto wrap = (std::forward<F>(f)(), detail::wrap_void());
            pwrap = wrap_into_unique(std::move(wrap));
        }
    }

    return pwrap->get();
}

template <typename T>
class FastLocalTransferWrap {
public:
    FastLocalTransferWrap(T &v, int root = 0) : p_{&v}, root_{root} {}

    T &ref() {
        return *p_;
    }

    const T &ref() const {
        return *p_;
    }

    void BinWrite(std::ostream &os) const {
        if (partask::master(root_)) {
            T *p = new T(std::move(const_cast<T&>(ref())));
            DEBUG("Writing address " << p);
            os.put('\1');
            io::binary::BinWrite(os, p);
        } else {
            os.put('\0');
            io::binary::BinWrite(os, ref());
        }
    }

    void BinRead(std::istream &is) {
        if (is.get()) {
            VERIFY(partask::master(root_));
            T *p = io::binary::BinRead<T*>(is);
            DEBUG("Readind address " << p);
            ref() = std::move(*p);
            delete p;
        } else {
            io::binary::BinRead(is, ref());
        }
    }

private:
    T *p_;
    int root_;
};

template <typename T>
auto fast_local_transfer(T &v, int root = 0) {
    return FastLocalTransferWrap<T>(v, root);
}

inline auto chunks_rr(size_t sz) {
    std::vector<size_t> chunks;
    size_t mpi_size = world_size();
    size_t mpi_rank = world_rank();
    for (size_t i = 0; i < sz; ++i) {
        if (i % mpi_size == mpi_rank)
            chunks.push_back(i);
    }

    return chunks;
}

template<class StreamListType>
void swap_streams(StreamListType &all_streams,
                  StreamListType &streams,
                  const std::vector<size_t> &chunks) {
    VERIFY(streams.size() == chunks.size());
    for (size_t i = 0; i < chunks.size(); ++i) {
        DEBUG("Swapping: " << i << " <-> " << chunks[i]);
        std::swap(streams[i], all_streams[chunks[i]]);
    }
}

template<class StreamListType>
auto create_empty_stream_list(size_t size) {
    StreamListType streams;
    for (size_t i = 0; i < size; ++i) {
        streams.push_back({});
    }
    return streams;
}

template<class StreamListType, class F>
void execute_on_subset(StreamListType &all_streams,
                       F f) {
    // Select streams
    std::vector<size_t> chunks = partask::chunks_rr(all_streams.size());
    INFO("Selected streams: " << chunks);
    execute_on_subset(all_streams, chunks, std::move(f));
}

template<class StreamListType, class F>
void execute_on_subset(StreamListType &all_streams,
                       const std::vector<size_t> &chunks,
                       F f) {
    auto local_streams = partask::create_empty_stream_list<StreamListType>(chunks.size());
    partask::swap_streams(all_streams, local_streams, chunks);
    local_streams.reset();
    f(local_streams);
    local_streams.close();
    partask::swap_streams(all_streams, local_streams, chunks);
}

}  // namespace partask

namespace io {

namespace binary {

// Enables io::binary::BinRead(is, FastLocalTransferWrap(x));
// It is not a specialization!
template <typename T>
void BinRead(std::istream &is, partask::FastLocalTransferWrap<T> &&w) {
    w.BinRead(is);
}

}  // namespace binary
}  // namespace io
