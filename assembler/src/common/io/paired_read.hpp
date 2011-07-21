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
 * PairedRead is a structure, where information from input files is stored.
 * It includes 2 SingleRead elements and a distance between them.
 */

#ifndef COMMON_IO_PAIREDREAD_HPP_
#define COMMON_IO_PAIREDREAD_HPP_

#include <string>
#include <utility>
#include "common/io/single_read.hpp"

class PairedRead {
 public:
  /*
   * Type of variables which will store file names for reading from
   * Reader stream.
   */ 
  typedef std::pair<std::string, std::string> FilenameType;

  /*
   * Default constructor.
   */
  PairedRead() {}

  /*
   * Conctructor from single reads.
   *
   * @param first First single read in the pair.
   * @param second Second single read in the pair.
   * @param distance Distance between two single reads.
   */
  PairedRead(const SingleRead& first,
             const SingleRead& second,
             size_t distance)
    : first_(first), second_(second), distance_(distance) {}

  /*
   * Return first single read in the pair.
   *
   * @return First single read.
   */
  const SingleRead& first() const {
    return first_;
  }

  /*
   * Return second single read in the pair.
   *
   * @return Second single read.
   */
  const SingleRead& second() const {
    return second_;
  }

  /*
   * Return distance of paired read.
   *
   * @return Distance.
   */
  size_t distance() const {
    return distance_;
  }

  /*
   * Check whether paired read is valid.
   *
   * @return true if paired read is valid (both single reads are correct),
   * and false otherwise.
   */
  bool IsValid() const {
    return first_.IsValid() && second_.IsValid();
  }

  /*
   * Return ith single read of pair (0th or 1st). If index
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
    assert(false);
  }

  /*
   * Return reversed complimentary paired read (paired read with
   * reserve complimentary first and second single reads
   * and the same distance.
   *
   * @return Reversed complimentary paired read.
   */
  const PairedRead operator!() const {
    return PairedRead(!second_, !first_, distance_);
  }

  /*
   * Check whether two paired reads are equal.
   *
   * @param pairedread The paired read we want to compare ours with.
   * @return true if these two paired reads have similar
   * first and second single reads and distance,
   * and false otherwise.
   */
  bool operator==(const PairedRead& pairedread) const {
    return first_ == pairedread.first_ &&
      second_ == pairedread.second_ &&
      distance_ == pairedread.distance_;
  }

 private:
  /*
   * @variable First single read in the pair.
   */
  SingleRead first_;
  /*
   * @variable Second single read in the pair.
   */
  SingleRead second_;
  /*
   * @variable Distance between two single reads.
   */
  size_t distance_;
};

#endif /* COMMON_IO_PAIREDREAD_HPP_ */
