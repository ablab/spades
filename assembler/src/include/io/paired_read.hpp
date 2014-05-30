//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "single_read.hpp"

#include <string>
#include <utility>

namespace io {

/**
 * It includes 2 SingleRead elements and the insert size.
 */
class PairedRead {
 public:
  typedef SingleRead SingleReadT;
  typedef int16_t size_type;

  /*
   * Default constructor.
   */
  PairedRead() : first_(), second_(), insert_size_(0) {}

  /*
   * Conctructor from SingleReads.
   *
   * @param first First SingleRead in the pair.
   * @param second Second SingleRead in the pair.
   * @param insert_size Insert size of the paired read.
   */
  PairedRead(const SingleRead& first,
             const SingleRead& second,
             size_t insert_size)
      : first_(first), second_(second), insert_size_(insert_size) {}

  /*
   * Return first SingleRead in the pair.
   *
   * @return First SingleRead.
   */
  const SingleRead& first() const {
    return first_;
  }

  /*
   * Return second SingleRead in the pair.
   *
   * @return Second SingleRead.
   */
  const SingleRead& second() const {
    return second_;
  }

  /*
   * Return insert_size of PairedRead.
   *
   * @return Insert size.
   */
  size_t insert_size() const {
    return insert_size_;
  }

  /*
   * Return distance between single reads.
   *
   * @return Distance.
   */
  size_t distance() const {
    return insert_size_ - second_.size();
  }

  /*
   * Return gap between single reads.
   *
   * @return Gap.
   */
  size_t gap() const {
    return insert_size_ - first_.size() - second_.size();
  }

  size_t size() const {
    return std::max(first_.size(), second_.size());
  }

  size_t nucl_count() const {
    return first_.size() + second_.size();
  }

  /*
   * Check whether PairedRead is valid.
   *
   * @return true if PairedRead is valid (both SingleReads are
   * correct), and false otherwise.
   */
  bool IsValid() const {
    return first_.IsValid() && second_.IsValid();
  }

  /*
   * Return ith SingleRead of pair (0th or 1st). If index
   * is not 0 or 1, the assertion happens.
   *
   * @param i SingleRead index.
   * @return SingleRead on ith position of pair.
   */
  const SingleRead& operator[] (size_t i) const {
    if (i == 0) {
      return first_;
    } else if (i == 1) {
      return second_;
    }
    VERIFY(false);
    return first_;
  }

  /*
   * Return reversed complimentary PairedRead (PairedRead with
   * reserve complimentary first and second SingleReads
   * and the same insert size).
   *
   * @return Reversed complimentary PairedRead.
   */
  const PairedRead operator!() const {
    return PairedRead(!second_, !first_, insert_size_);
  }

  /*
   * Check whether two PairedReads are equal.
   *
   * @param pairedread The PairedRead we want to compare ours with.
   * @return true if these two PairedReads have similar
   * first and second SingleReads and insert size,
   * and false otherwise.
   */
  bool operator==(const PairedRead& pairedread) const {
    return first_ == pairedread.first_ &&
        second_ == pairedread.second_ &&
        insert_size_ == pairedread.insert_size_;
  }

  bool BinWrite(std::ostream& file, bool rc1 = false, bool rc2 = false) const {
    first_.BinWrite(file, rc1);
    second_.BinWrite(file, rc2);

    size_type is = (size_type)insert_size_;
    file.write((const char *) &is, sizeof(is));

    return !file.fail();
  }

  void print_size() const {
    first_.print_size();
    second_.print_size();
  }

 private:
  /*
   * @variable First SingleRead in the pair.
   */
  SingleRead first_;
  /*
   * @variable Second SingleRead in the pair.
   */
  SingleRead second_;
  /*
   * @variable Insert size between two SingleReads.
   */
  size_t insert_size_;

};

inline std::ostream& operator<<(std::ostream& os, const PairedRead& read) {
  os << "Single read first=" << read.first() << " second=" << read.second() << std::endl;
  return os;
}

class PairedReadSeq {
 public:
  typedef SingleReadSeq SingleReadT;
 private:
  SingleReadSeq first_;
  SingleReadSeq second_;
  size_t insert_size_;

 public:
  PairedReadSeq() : first_(), second_(), insert_size_(0) {}

  //    PairedReadSeq(std::istream& file, size_t is): first_(file), second_(file) {
  //        PairedRead::size_type is_delta;
  //        file.read((char *) &is_delta, sizeof(is_delta));
  //
  //        insert_size_ = is + is_delta;
  //    }

  bool BinRead(std::istream& file, size_t is = 0) {
    first_.BinRead(file);
    second_.BinRead(file);

    PairedRead::size_type is_delta;
    file.read((char *) &is_delta, sizeof(is_delta));

    insert_size_ = is + is_delta;
    return !file.fail();
  }

  bool BinWrite(std::ostream& file) const {
    first_.BinWrite(file);
    second_.BinWrite(file);

    PairedRead::size_type is = (PairedRead::size_type)insert_size_;
    file.write((const char *) &is, sizeof(is));

    return !file.fail();
  }

  const SingleReadSeq& first() const {
    return first_;
  }

  const SingleReadSeq& second() const {
    return second_;
  }

  size_t insert_size() const {
    return insert_size_;
  }

  size_t distance() const {
    return insert_size_ - second_.size();
  }

  size_t gap() const {
    return insert_size_ - first_.size() - second_.size();
  }

  size_t size() const {
    return std::max(first_.size(), second_.size());
  }

  size_t nucl_count() const {
    return first_.size() + second_.size();
  }

  PairedReadSeq(const SingleReadSeq& first,
                const SingleReadSeq& second,
                size_t insert_size)
      : first_(first), second_(second), insert_size_(insert_size) {}

  const SingleReadSeq& operator[] (size_t i) const {
    if (i == 0) {
      return first_;
    } else if (i == 1) {
      return second_;
    }
    VERIFY(false);
    return first_;
  }

  const PairedReadSeq operator!() const {
    return PairedReadSeq(!second_, !first_, insert_size_);
  }

};

inline std::ostream& operator<<(std::ostream& os, const PairedReadSeq& read) {
  os << "Paired read first=" << read.first() << " second=" << read.second() << std::endl;
  return os;
}

}
