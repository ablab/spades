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

#ifndef PAIREDREAD_HPP_
#define PAIREDREAD_HPP_

#include "common/io/single_read.hpp"

class PairedRead {
private:
	SingleRead first_;
	SingleRead second_;
	size_t distance_;
public:

	PairedRead() {}

	PairedRead(const SingleRead& first, const SingleRead& second, size_t distance) 
    : first_(first), second_(second), distance_(distance) {}

	const SingleRead& first() const {
		return first_;
	}

	const SingleRead& second() const {
		return second_;
	}

	size_t distance() const {
		return distance_;
	}

	bool IsValid() const {
		return first_.IsValid() && second_.IsValid();
	}

	const SingleRead& operator[] (size_t index) const {
		if (index == 0) {
			return first_;
		} else if (index == 1) {
			return second_;
		}
		assert(false);
	}
  
	const PairedRead operator!() const {
		return PairedRead(!second_, !first_, distance_);
  }

  bool operator==(const PairedRead& pairedread) const {
    return first_ == pairedread.first_ &&
      second_ == pairedread.second_ &&
      distance_ == pairedread.distance_;
  }
};

#endif /* PAIRED_READ_HPP_ */
