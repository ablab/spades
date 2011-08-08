/*
 * math.h
 *
 *  Created on: Jul 30, 2011
 *      Author: sergey
 */

#ifndef XMATH_H_
#define XMATH_H_

#include "sequence/seq.hpp"

namespace math
{
    template<class T>
    T eps();

    template<>
    inline double eps<double>() { return 1e-10; }

    template<>
    inline float eps<float>() { return 1e-5; }

	template<class T>
	bool eq(T lhs, T rhs) { return abs(lhs - rhs) < eps<T>(); }

	template<class T>
	bool ls(T lhs, T rhs) { return lhs + eps<T>() < rhs; }

	template<class T>
	bool gr(T lhs, T rhs) { return ls(rhs, lhs); }

	template<class T>
	bool le(T lhs, T rhs) { return !gr(lhs, rhs); }

	template<class T>
	bool ge(T lhs, T rhs) { return !ls(lhs, rhs); }
}


#endif /* XMATH_H_ */
