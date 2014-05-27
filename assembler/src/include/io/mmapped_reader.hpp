//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef HAMMER_MMAPPED_READER_HPP
#define HAMMER_MMAPPED_READER_HPP

#include "adt/pointer_iterator.hpp"
#include "adt/array_vector.hpp"

#include "verify.hpp"

#include <boost/iterator/iterator_facade.hpp>

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
  MMappedReader()
      : StreamFile(-1), Unlink(false), FileName(""), MappedRegion(0), FileSize(0), BytesRead(0)
    {}

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

  MMappedReader(MMappedReader &&other) {
    // First, copy out the stuff
    MappedRegion = other.MappedRegion;
    FileSize = other.FileSize;
    BlockOffset = other.BlockOffset;
    BytesRead = other.BytesRead;
    BlockSize = other.BlockSize;
    FileName = std::move(other.FileName);
    Unlink = other.Unlink;
    StreamFile = other.StreamFile;

    // Now, zero out inside other, so we won't do crazy thing in dtor
    StreamFile = -1;
    Unlink = false;
    MappedRegion = 0;
  }

  MMappedReader& operator=(MMappedReader &&other) {
    if (this != &other) {
      *this = std::move(other);
    }
    return *this;
  }

  virtual ~MMappedReader() {
    if (StreamFile != -1)
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

  void* skip(size_t amount) {
    // Easy case, no remapping is needed
    if (BytesRead + amount <= BlockOffset + BlockSize) {
      void* out = MappedRegion + BytesRead - BlockOffset;
      BytesRead += amount;

      return out;
    }

    // Make sure data does not cross the block boundary
    VERIFY(BytesRead  == BlockOffset + BlockSize);

    // Now, remap and read from the beginning of the block
    remap();

    return skip(amount);
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
                      size_t blocksize = 64*1024*1024 / (sizeof(T) * (unsigned)getpagesize()) * (sizeof(T) * (unsigned)getpagesize()),
                      size_t off = 0, size_t sz = 0):
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

template<class T>
class MMappedFileRecordIterator :
        public boost::iterator_facade<MMappedFileRecordIterator<T>,
                                      const T,
                                      std::input_iterator_tag> {
  public:
    // Default ctor, used to implement "end" iterator
    MMappedFileRecordIterator() : good_(false) {}
    MMappedFileRecordIterator(const std::string &FileName)
            : reader_(FileName, false), good_(true) {
        reader_.read(&value_, sizeof(value_));
    }
    MMappedFileRecordIterator(MMappedRecordReader<T> &&reader)
            : reader_(std::move(reader)), good_(true) {
        reader_.read(&value_, sizeof(value_));
    }
    bool good() const {
        return good_;
    }

  private:
    friend class boost::iterator_core_access;

    void increment() {
        good_ = reader_.good();
        if (good_)
            reader_.read(&value_, sizeof(value_));
    }
    bool equal(const MMappedFileRecordIterator &other) {
        // Iterators are equal iff:
        //   1) They both are not good (at the end of the stream),
        //      or
        //   2) Has the same mapped region
        return ((!reader_.good() && !other.reader_.good()) ||
                reader_.data() == other.reader_.data());
    }
    const T dereference() const { return value_; }

    T value_;
    MMappedRecordReader<T> reader_;
    bool good_;
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

template<class T>
class MMappedFileRecordArrayIterator :
        public boost::iterator_facade<MMappedFileRecordArrayIterator<T>,
                                      const T*,
                                      std::input_iterator_tag,
                                      const T*> {
  public:
    // Default ctor, used to implement "end" iterator
    MMappedFileRecordArrayIterator(): value_(NULL), reader_(), elcnt_(0), good_(false) {}
    MMappedFileRecordArrayIterator(const std::string &FileName, size_t elcnt)
            : value_(NULL),
              reader_(FileName, false,
                      64*1024*1024 / (sizeof(T) * (unsigned)getpagesize() * elcnt) * (sizeof(T) * (unsigned)getpagesize() * elcnt)),
              elcnt_(elcnt), good_(true) {
        increment();
    }
    MMappedFileRecordArrayIterator(MMappedRecordReader<T> &&reader, size_t elcnt)
            : value_(NULL), reader_(std::move(reader)), elcnt_(elcnt), good_(true) {
        increment();
    }

    bool good() const { return good_; }

  private:
    friend class boost::iterator_core_access;

    void increment() {
        good_ = reader_.good();
        if (good_)
            value_ = (T*)reader_.skip(elcnt_ * sizeof(T));
    }
    bool equal(const MMappedFileRecordArrayIterator &other) const {
        return value_ == other.value_;
    }
    const T* dereference() const { return value_; }

    T* value_;
    MMappedRecordReader<T> reader_;
    size_t elcnt_;
    bool good_;
};

#endif // HAMMER_MMAPPED_READER_HPP
