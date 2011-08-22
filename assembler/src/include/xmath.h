/*
 * xmath.h
 *
 *  Created on: Jul 30, 2011
 *      Author: sergey
 */

#ifndef XMATH_H_
#define XMATH_H_

#include <cmath>

namespace math
{
    template<class T>
    T eps();

    template<>
    inline double eps<double>() { return 1e-10; }

    template<>
    inline float eps<float>() { return 1e-5; }

	template<class T>
	bool ls(T lhs, T rhs) { return lhs + eps<T>() < rhs; }

	//todo return initial variant fixing compilation
	template<class T>
	bool eq(T lhs, T rhs) { return !ls(lhs, rhs) && !ls(rhs, lhs)/* std::abs(lhs - rhs) < eps<T>()*/; }

	template<class T>
	bool gr(T lhs, T rhs) { return ls(rhs, lhs); }

	template<class T>
	bool le(T lhs, T rhs) { return !gr(lhs, rhs); }

	template<class T>
	bool ge(T lhs, T rhs) { return !ls(lhs, rhs); }

	template<class T>
	T floor(T t) {return std::floor(t + eps<T>());}

	template<class T>
	T round(T t) {return floor(t + 0.5);}

}


#endif /* XMATH_H_ */
