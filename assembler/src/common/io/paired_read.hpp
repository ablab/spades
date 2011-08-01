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
  PairedRead() : first_(), second_(), distance_(0) {}

  /*
   * Conctructor from SingleReads.
   *
   * @param first First SingleRead in the pair.
   * @param second Second SingleRead in the pair.
   * @param distance Distance between two SingleReads.
   */
  PairedRead(const SingleRead& first,
             const SingleRead& second,
             size_t distance)
    : first_(first), second_(second), distance_(distance) {}

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
   * Return distance of PairedRead.
   *
   * @return Distance.
   */
  size_t distance() const {
    return distance_;
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
    assert(false);
  }

  /*
   * Return reversed complimentary PairedRead (PairedRead with
   * reserve complimentary first and second SingleReads
   * and the same distance).
   *
   * @return Reversed complimentary PairedRead.
   */
  const PairedRead operator!() const {
    return PairedRead(!second_, !first_, distance_);
  }

  /*
   * Check whether two PairedReads are equal.
   *
   * @param pairedread The PairedRead we want to compare ours with.
   * @return true if these two PairedReads have similar
   * first and second SingleReads and distance,
   * and false otherwise.
   */
  bool operator==(const PairedRead& pairedread) const {
    return first_ == pairedread.first_ &&
      second_ == pairedread.second_ &&
      distance_ == pairedread.distance_;
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
   * @variable Distance between two SingleReads.
   */
  size_t distance_;
};

#endif /* COMMON_IO_PAIREDREAD_HPP_ */
