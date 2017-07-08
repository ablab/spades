//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef HAMMER_MMAPPED_WRITER_HPP
#define HAMMER_MMAPPED_WRITER_HPP

#include "adt/pointer_iterator.hpp"
#include "adt/array_vector.hpp"
#include "common/utils/verify.hpp"

#include <string>

#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <strings.h>
#include <stdio.h>

class MMappedWriter {
    int StreamFile;

    MMappedWriter(const MMappedWriter &) = delete;

protected:
    uint8_t *MappedRegion;
    size_t BytesWritten, BytesReserved, FileOffset, BufOffset;
public:
    MMappedWriter() = default;

    MMappedWriter(const std::string &FileName) {
        open(FileName);
    }

    void open(const std::string &FileName) {
        StreamFile = ::open(FileName.c_str(), O_RDWR | O_CREAT | O_TRUNC, (mode_t) 0660);
        VERIFY_MSG(StreamFile != -1,
                   "open(2) failed. Reason: " << strerror(errno) << ". Error code: " << errno);

        FileOffset = BytesWritten = 0;
        MappedRegion = NULL;
    }

    virtual ~MMappedWriter() {
        if (MappedRegion)
            munmap(MappedRegion, BytesReserved);
        close(StreamFile);
    }

    void write(void *buf, size_t amount) {
        memcpy(MappedRegion + BufOffset + BytesWritten, buf, amount);
        BytesWritten += amount;
    }

    bool good() const {
        return BytesWritten < BytesReserved;
    }

    void reserve(size_t amount) {
        if (MappedRegion) {
            munmap(MappedRegion, BytesReserved);
            FileOffset += BytesWritten;
            MappedRegion = NULL;
        }

        if (amount == 0)
            return;

        int res = (int) lseek(StreamFile, amount - 1, SEEK_CUR);
        VERIFY_MSG(res != -1,
                   "lseek(2) failed. Reason: " << strerror(errno) << ". Error code: " << errno);
        res = (int) ::write(StreamFile, "", 1);
        VERIFY_MSG(res != -1,
                   "write(2) failed. Reason: " << strerror(errno) << ". Error code: " << errno);

        // FileOffset here should be aligned to page boundary. Tune the stuff due to this fact.
        int PageSize = getpagesize();
        size_t FileOffsetAligned = FileOffset / PageSize * PageSize;
        size_t Residual = FileOffset - FileOffsetAligned;

        BytesReserved = amount + Residual;
        BytesWritten = 0;
        BufOffset = Residual;
        MappedRegion =
                (uint8_t *) mmap(NULL, BytesReserved,
                                 PROT_READ | PROT_WRITE, MAP_FILE | MAP_SHARED,
                                 StreamFile, FileOffsetAligned);
        VERIFY_MSG((intptr_t) MappedRegion != -1L,
                   "mmap(2) failed. Reason: " << strerror(errno) << ". Error code: " << errno);
    }

    size_t size() const { return BytesReserved; }
};

template<typename T>
class MMappedRecordWriter : public MMappedWriter {
public:
    typedef adt::pointer_iterator<T> iterator;
    typedef const adt::pointer_iterator<T> const_iterator;

    MMappedRecordWriter() = default;

    MMappedRecordWriter(const std::string &FileName) :
            MMappedWriter(FileName) {
    }

    void write(const T *el, size_t amount) {
        MMappedWriter::write((void *) el, amount * sizeof(T));
    }

    void reserve(size_t amount) {
        MMappedWriter::reserve(amount * sizeof(T));
    }

    void resize(size_t amount) {
        MMappedWriter::reserve(amount * sizeof(T));
    }

    size_t size() const { return BytesReserved / sizeof(T); }

    T *data() { return (T *) MappedRegion; }

    const T *data() const { return (const T *) MappedRegion; }

    T &operator[](size_t idx) { return data()[idx]; }

    const T &operator[](size_t idx) const { return data()[idx]; }

    iterator begin() { return iterator(data()); }

    const_iterator begin() const { return const_iterator(data()); }

    iterator end() { return iterator(data() + size()); }

    const_iterator end() const { return const_iterator(data() + size()); }
};

template<typename T>
class MMappedRecordArrayWriter : public MMappedWriter {
    size_t elcnt_;
public:
    typedef typename adt::array_vector<T>::iterator iterator;
    typedef typename adt::array_vector<T>::const_iterator const_iterator;

    MMappedRecordArrayWriter() = default;

    MMappedRecordArrayWriter(const std::string &FileName,
                             size_t elcnt = 1) :
            MMappedWriter(FileName), elcnt_(elcnt) { }

    void open(const std::string &FileName,
              size_t elcnt = 1) {
        elcnt_ = elcnt;
        MMappedWriter::open(FileName);
    }

    void write(const T *el, size_t amount) {
        MMappedWriter::write((void *) el, amount * sizeof(T) * elcnt_);
    }

    void reserve(size_t amount) {
        MMappedWriter::reserve(amount * sizeof(T) * elcnt_);
    }

    void resize(size_t amount) {
        MMappedWriter::reserve(amount * sizeof(T) * elcnt_);
    }

    size_t size() const { return BytesReserved / sizeof(T) / elcnt_; }

    T *data() { return (T *) MappedRegion; }

    const T *data() const { return (const T *) MappedRegion; }

    T &operator[](size_t idx) { return data()[idx * elcnt_]; }

    const T &operator[](size_t idx) const { return data()[idx * elcnt_]; }

    iterator begin() { return iterator(data(), elcnt_); }

    const_iterator begin() const { return const_iterator(data(), elcnt_); }

    iterator end() { return iterator(data() + size() * elcnt_, elcnt_); }

    const_iterator end() const { return const_iterator(data() + size() * elcnt_, elcnt_); }
};

#endif // HAMMER_MMAPPED_WRITER_HPP
