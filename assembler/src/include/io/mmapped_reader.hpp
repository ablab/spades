//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef HAMMER_MMAPPED_READER_HPP
#define HAMMER_MMAPPED_READER_HPP

#include "pointer_iterator.hpp"

#include "verify.hpp"

#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <strings.h>
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
    VERIFY((intptr_t)MappedRegion != -1L);
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
                size_t blocksize = 64*1024*1024)
      : Unlink(unlink), FileName(filename), BlockSize(blocksize) {
    struct stat buf;

    FileSize = (stat(FileName.c_str(), &buf) != 0 ? 0 : buf.st_size);

    StreamFile = open(FileName.c_str(), O_RDONLY);
    VERIFY(StreamFile != -1);

    if (BlockSize != -1ULL) {
      size_t PageSize = getpagesize();
      BlockSize = BlockSize / PageSize * PageSize;
    } else
      BlockSize = FileSize;

    if (BlockSize) {
      MappedRegion =
          (uint8_t*)mmap(NULL, BlockSize, PROT_READ | PROT_WRITE, MAP_FILE | MAP_PRIVATE,
                         StreamFile, 0);
      VERIFY((intptr_t)MappedRegion != -1L);
    } else
      MappedRegion = NULL;

    BlockOffset = BytesRead = 0;
  }

  virtual ~MMappedReader() {
    if (MappedRegion)
      munmap(MappedRegion, BlockSize);
    close(StreamFile);

    if (Unlink) {
      int res = unlink(FileName.c_str());
      VERIFY(res == 0);
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

  void* data() const { return MappedRegion; }
};

template<typename T>
class MMappedRecordReader : public MMappedReader {
 public:
  typedef pointer_iterator<T> iterator;
  typedef const pointer_iterator<T> const_iterator;
  
  MMappedRecordReader(const std::string &FileName, bool unlink = true,
                      size_t blocksize = 64*1024*1024):
      MMappedReader(FileName, unlink, blocksize) {
    VERIFY(FileSize % sizeof(T) == 0);
  }

  void read(T* el, size_t amount) {
    MMappedReader::read(el, amount * sizeof(T));
  }

  size_t size() const { return FileSize / sizeof(T); }
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
  typedef pointer_array_iterator<T> iterator;
  typedef const pointer_array_iterator<T> const_iterator;

  MMappedRecordArrayReader(const std::string &FileName,
                           size_t elcnt = 1,
                           bool unlink = true,
                           size_t blocksize = 64*1024*1024):
      MMappedReader(FileName, unlink, blocksize), elcnt_(elcnt){
    VERIFY(FileSize % (sizeof(T) * elcnt_) == 0);
  }

  void read(T* el, size_t amount) {
    MMappedReader::read(el, amount * sizeof(T) * elcnt_);
  }

  size_t size() const { return FileSize / sizeof(T) / elcnt_; }
  T* data() { return (T*)MappedRegion; }
  const T* data() const { return (const T*)MappedRegion; }
  T& operator[](size_t idx) { return data()[idx*elcnt_]; }
  const T& operator[](size_t idx) const { return data()[idx*elcnt_]; }

  iterator begin() { return iterator(data(), /* size */ elcnt_); }
  const_iterator begin() const { return const_iterator(data()), /* size */ elcnt_; }
  iterator end() { return iterator(data() + size()*elcnt_, elcnt_); }
  const_iterator end() const { return const_iterator(data() + size()*elcnt_, elcnt_); }
};


#endif // HAMMER_MMAPPED_READER_HPP
