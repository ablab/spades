/*
 * normalize.hpp
 *
 *  Created on: Sep 19, 2011
 *      Author: andrey
 */

#ifndef NORMALIZE_HPP_
#define NORMALIZE_HPP_

#include "lc_common.hpp"

namespace long_contigs {

class Normalizer {
public:
	Normalizer(size_t insert_szie, size_t read_length): is(insert_szie), rs(read_length) {
	}

	//Weight of ideal info
	virtual int IdealWeight(Graph& g, EdgeId edge1, EdgeId edge2, int distance, size_t delta = 0) = 0;

	virtual double Normalize(double weight, Graph& g, EdgeId edge1, EdgeId edge2, size_t distance, size_t delta = 0) = 0;

protected:
	size_t is;
	size_t rs;
};


class ReadCountNormalizer: public Normalizer {
public:

	virtual int IdealWeight(Graph& g, EdgeId edge1, EdgeId edge2, int distance, size_t delta = 0) {
		if (distance < 0) {
			return IdealWeigth(g, e2, e1, -distance, delta);
		}
		int gaplLen = distance - g.length(e1);

		int right = std::min(is, gaplLen + g.length(e2) + rs);
		int left = std::max(gaplLen - is, - rs - g.length(e1)) + is;

		int delta = right - left + 1 - debruijn::K + delta;

		return delta;
	}
};

} //namespace long_contigs


#endif /* NORMALIZE_HPP_ */
