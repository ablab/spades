//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef HAMMER_MMAPPED_WRITER_HPP
#define HAMMER_MMAPPED_WRITER_HPP

#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <strings.h>
#include <string>

class MMappedWriter {
  int StreamFile;
  MMappedWriter(const MMappedWriter &) = delete;
 protected:
  uint8_t* MappedRegion;
  size_t BytesWritten, BytesReserved, FileOffset, BufOffset;
 public:
  MMappedWriter() = default;
  MMappedWriter(const std::string &FileName) {
    open(FileName);
  }

  void open(const std::string &FileName) {
    StreamFile = ::open(FileName.c_str(), O_RDWR | O_CREAT | O_TRUNC, (mode_t)0660);
    VERIFY(StreamFile != -1);
    
    FileOffset = BytesWritten = 0;
    MappedRegion = NULL;
  }
  
  virtual ~MMappedWriter() {
    if (MappedRegion)
      munmap(MappedRegion, BytesReserved);
    close(StreamFile);
  }

  void write(void* buf, size_t amount) {
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
    }
    
    int res = lseek(StreamFile, amount-1, SEEK_CUR);
    VERIFY(res != -1);
    res = ::write(StreamFile, "", 1);
    VERIFY(res != -1);
    
    // FileOffset here should be aligned to page boundary. Tune the stuff due to this fact.
    int PageSize = getpagesize();
    size_t FileOffsetAligned = FileOffset / PageSize * PageSize;
    size_t Residual = FileOffset - FileOffsetAligned;
    
    BytesReserved = amount + Residual;
    BytesWritten = 0; BufOffset = Residual;
    MappedRegion =
        (uint8_t*)mmap(NULL, BytesReserved,
                       PROT_READ | PROT_WRITE, MAP_FILE | MAP_SHARED,
                       StreamFile, FileOffsetAligned);
    VERIFY((intptr_t)MappedRegion != -1L);
  }
  
  size_t size() const { return BytesReserved; }
};

template<typename T>
class MMappedRecordWriter : public MMappedWriter {
 public:
  typedef pointer_iterator<T> iterator;
  typedef const pointer_iterator<T> const_iterator;

  MMappedRecordWriter() = default;
  MMappedRecordWriter(const std::string &FileName):
      MMappedWriter(FileName) {
  }

  void write(const T* el, size_t amount) {
    MMappedWriter::write((void*)el, amount * sizeof(T));
  }

  void reserve(size_t amount) {
    MMappedWriter::reserve(amount * sizeof(T));
  }

  void resize(size_t amount) {
    MMappedWriter::reserve(amount * sizeof(T));
  }

  size_t size() const { return BytesReserved / sizeof(T); }
  T* data() { return (T*)MappedRegion; }
  const T* data() const { return (const T*)MappedRegion; }
  T& operator[](size_t idx) { return data()[idx]; }
  const T& operator[](size_t idx) const { return data()[idx]; }

  iterator begin() { return iterator(data()); }
  const_iterator begin() const { return const_iterator(data()); }
  iterator end() { return iterator(data()+ size()); }
  const_iterator end() const { return const_iterator(data() + size()); }
};

#endif // HAMMER_MMAPPED_WRITER_HPP
