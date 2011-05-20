/*
 * paired_read.hpp
 *
 *  Created on: 19.05.2011
 *      Author: vyahhi
 */

#ifndef PAIRED_READ_HPP_
#define PAIRED_READ_HPP_

class PairedRead {
private:
	Read first_;
	Read second_;
	size_t distance_;
public:

	PairedRead()  {

	}

	PairedRead(const Read& first, const Read& second, size_t distance) : first_(first), second_(second), distance_(distance) {

	}

	const Read& first() const {
		return first_;
	}

	const Read& second() const {
		return second_;
	}

	size_t distance() const {
		return distance_;
	}

	bool isValid() const {
		return first_.isValid() && second_.isValid();
	}

	const Read& operator[] (size_t index) const {
		if (index == 0) {
			return first_;
		} else if (index == 1) {
			return second_;
		}
		assert(false);
	}

	const PairedRead operator!() const{
		return PairedRead(!second_, !first_, distance_);
	}

};

// todo: put this to *.cpp
//ostream& operator<<(ostream& os, const PairedRead& p_r) {
//	return os << "[" << p_r.first() << ", " << p_r.second() << "]";
//}

#endif /* PAIRED_READ_HPP_ */
