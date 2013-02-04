//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef HAMMER_MMAPPED_READER_HPP
#define HAMMER_MMAPPED_READER_HPP

#include "adt/pointer_iterator.hpp"
#include "adt/array_vector.hpp"

#include "verify.hpp"

#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <cstring>
#include <cerrno>

#include <string>

class MMappedReader {
  int StreamFile;
  bool Unlink;
  std::string FileName;
  MMappedReader(const MMappedReader &) = delete;

  void remap() {
    VERIFY(BlockSize != FileSize);

    if (MappedRegion)
      munmap(MappedRegion, BlockSize);

    BlockOffset += BlockSize;
    // We do not add PROT_WRITE here intentionaly - remapping and write access
    // is pretty error-prone.
    MappedRegion =
        (uint8_t*)mmap(NULL, BlockSize,
                       PROT_READ, MAP_FILE | MAP_PRIVATE,
                       StreamFile, BlockOffset);
    VERIFY_MSG((intptr_t)MappedRegion != -1L,
               "mmap(2) failed. Reason: " << strerror(errno) << ". Error code: " << errno);
  }

  void read_internal(void *buf, size_t amount) {
    memcpy(buf, MappedRegion + BytesRead - BlockOffset, amount);
    BytesRead += amount;
  }


 protected:
  uint8_t* MappedRegion;
  size_t FileSize, BlockOffset, BytesRead, BlockSize;

 public:
  MMappedReader(const std::string &filename, bool unlink = false,
                size_t blocksize = 64*1024*1024, size_t off = 0, size_t sz = 0)
      : Unlink(unlink), FileName(filename), BlockSize(blocksize) {
    struct stat buf;

    FileSize = (sz ? sz : (stat(FileName.c_str(), &buf) != 0 ? 0 : buf.st_size));

    StreamFile = open(FileName.c_str(), O_RDONLY);
    VERIFY_MSG(StreamFile != -1,
               "open(2) failed. Reason: " << strerror(errno) << ". Error code: " << errno);

    if (BlockSize != -1ULL) {
      size_t PageSize = getpagesize();
      BlockSize = BlockSize / PageSize * PageSize;
    } else
      BlockSize = FileSize;

    if (BlockSize) {
      MappedRegion =
          (uint8_t*)mmap(NULL, BlockSize, PROT_READ | PROT_WRITE, MAP_FILE | MAP_PRIVATE,
                         StreamFile, off);
      VERIFY_MSG((intptr_t)MappedRegion != -1L,
                 "mmap(2) failed. Reason: " << strerror(errno) << ". Error code: " << errno);
    } else
      MappedRegion = NULL;

    BlockOffset = BytesRead = 0;
  }

  virtual ~MMappedReader() {
    close(StreamFile);
    if (MappedRegion)
      munmap(MappedRegion, BlockSize);

    if (Unlink) {
      int res = unlink(FileName.c_str());
      VERIFY_MSG(res == 0,
                 "unlink(2) failed. Reason: " << strerror(errno) << ". Error code: " << errno);
    }
  }

  void read(void* buf, size_t amount) {
    if (BytesRead + amount < BlockOffset + BlockSize) {
      // Easy case, no remap is necessary
      read_internal(buf, amount);
      return;
    }

    // Hard case - remapping is necessary. First - finish the current block.
    size_t ToRead = BlockSize - (BytesRead - BlockOffset);
    uint8_t *cbuf = (uint8_t*)buf;

    read_internal(cbuf, ToRead);
    amount -= ToRead; cbuf += ToRead;

    // Next, read as much BlockSize blocks as possible.
    while (amount >= BlockSize) {
      remap();
      read_internal(cbuf, BlockSize);
      amount -= BlockSize; cbuf += BlockSize;
    }

    // Finally, remap and read remaining.
    remap();
    read_internal(cbuf, amount);
  }

  bool good() const {
    return BytesRead < FileSize;
  }

  size_t size() const { return FileSize; }
  size_t data_size() const { return FileSize; }

  void* data() const { return MappedRegion; }
};

template<typename T>
class MMappedRecordReader : public MMappedReader {
 public:
  typedef pointer_iterator<T> iterator;
  typedef const pointer_iterator<T> const_iterator;

  MMappedRecordReader(const std::string &FileName, bool unlink = true,
                      size_t blocksize = 64*1024*1024, size_t off = 0, size_t sz = 0):
      MMappedReader(FileName, unlink, blocksize, off, sz) {
    VERIFY(FileSize % sizeof(T) == 0);
  }

  void read(T* el, size_t amount) {
    MMappedReader::read(el, amount * sizeof(T));
  }

  size_t size() const { return FileSize / sizeof(T); }
  size_t data_size() const { return FileSize; }
  T* data() { return (T*)MappedRegion; }
  const T* data() const { return (const T*)MappedRegion; }
  T& operator[](size_t idx) { return data()[idx]; }
  const T& operator[](size_t idx) const { return data()[idx]; }

  iterator begin() { return iterator(data()); }
  const_iterator begin() const { return const_iterator(data()); }
  iterator end() { return iterator(data()+ size()); }
  const_iterator end() const { return const_iterator(data() + size()); }
};

template<typename T>
class MMappedRecordArrayReader : public MMappedReader {
  size_t elcnt_;

 public:
  typedef typename array_vector<T>::iterator iterator;
  typedef typename array_vector<T>::const_iterator const_iterator;

  MMappedRecordArrayReader(const std::string &FileName,
                           size_t elcnt = 1,
                           bool unlink = true,
                           size_t off = 0, size_t sz = 0):
      MMappedReader(FileName, unlink, -1ULL, off, sz), elcnt_(elcnt) {
    VERIFY(FileSize % (sizeof(T) * elcnt_) == 0);
  }

  void read(T* el, size_t amount) {
    MMappedReader::read(el, amount * sizeof(T) * elcnt_);
  }

  size_t size() const { return FileSize / sizeof(T) / elcnt_; }
  size_t data_size() const { return FileSize; }
  size_t elcnt() const { return elcnt_; }
  T* data() { return (T*)MappedRegion; }
  const T* data() const { return (const T*)MappedRegion; }
  T& operator[](size_t idx) { return data()[idx*elcnt_]; }
  const T& operator[](size_t idx) const { return data()[idx*elcnt_]; }

  iterator begin() { return iterator(data(), /* size */ elcnt_); }
  const_iterator begin() const { return const_iterator(data()), /* size */ elcnt_; }
  const_iterator cbegin() const { return const_iterator(data()), /* size */ elcnt_; }
  iterator end() { return iterator(data() + size()*elcnt_, elcnt_); }
  const_iterator end() const { return const_iterator(data() + size()*elcnt_, elcnt_); }
  const_iterator cend() const { return const_iterator(data() + size()*elcnt_, elcnt_); }
};


#endif // HAMMER_MMAPPED_READER_HPP
