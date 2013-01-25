//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/**
 * @file    paired_read.hpp
 * @author  Mariya Fomkina
 * @version 1.0
 *
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 *
 * @section DESCRIPTION
 *
 * PairedRead is a structure, where information from input files is
 * stored. 
 * It includes 2 SingleRead elements and the insert size.
 */

#ifndef COMMON_IO_PAIREDREAD_HPP_
#define COMMON_IO_PAIREDREAD_HPP_

#include <string>
#include <utility>
#include "io/single_read.hpp"

namespace io {

class PairedRead {
 public:
  /*
   * Type of variables which will store file names for reading from
   * Reader stream.
   */ 
  typedef std::pair<std::string, std::string> FilenamesType;


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
      return max(first_.size(), second_.size());
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

  bool BinWrite(std::ostream& file) const {
      first_.BinWrite(file);
      second_.BinWrite(file);


      size_type is = insert_size_;
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

inline ostream& operator<<(ostream& os, const PairedRead& read) {
	os << "Single read first=" << read.first() << " second=" << read.second() << endl;
	return os;
}

class PairedReadSeq {

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

        PairedRead::size_type is = insert_size_;
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
        return max(first_.size(), second_.size());
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

    void inc_insert_size(size_t val) {
        insert_size_ += val;
    }

};

inline ostream& operator<<(ostream& os, const PairedReadSeq& read) {
	os << "Paired read first=" << read.first() << " second=" << read.second() << endl;
	return os;
}

}

#endif /* COMMON_IO_PAIREDREAD_HPP_ */
