/*
 * mathfunctions.hpp
 *
 *  Created on: 11.05.2011
 *      Author: snikolenko
 */

#ifndef MATHFUNCTIONS_H
#define MATHFUNCTIONS_H

#include<cmath>
#include<math.h>
#include "hammerread.hpp"

/**
  * @return logarithm of {n choose k}
  */
double logChooseNK(int n, int k) {
	double res = 0;
	for (int i=1; i <= k; ++i) {
		res += log(n-k+i) - log(i);
	}
	return res;
}

/**
  * @return Beta(x,y)
  */
double lBeta(int x, int y) {
	return (lgamma(x) + lgamma(y) - lgamma(x+y));
}

/**
  * @return log({a_1+...+a_n \choose a_1, ..., a_n})
  */
double lMultinomial(const vector<HammerRead> & x) {
	double res = 0.0, sum = 0.0;
	for (size_t i=0; i<x.size(); ++i) {
		res += lgamma(x[i].count+1);
		sum += x[i].count;
	}
	return (lgamma(sum+1) - res);
}

/**
  * @return log({a_1+...+a_n \choose a_1, ..., a_n}) for reads corresponding to the mask
  */
double lMultinomialWithMask(const vector<HammerRead> & x, const vector<int> & mask, int maskval) {
	assert(x.size() == mask.size());
	double res = 0.0, sum = 0.0;
	for (size_t i=0; i<x.size(); ++i) {
		if (mask[i] != maskval) continue;
		res += lgamma(x[i].count+1);
		sum += x[i].count;
	}
	return (lgamma(sum+1) - res);
}

/**
  * @return log(Beta(a_1+1, ..., a_n+1))
  */
double lBetaPlusOne(const vector<int> & x) {
	double res = 0.0, sum = 0.0;
	for (size_t i=0; i<x.size(); ++i) {
		res += lgamma(x[i]+1);
		sum += x[i]+1;
	}
	return (res - lgamma(sum));
}

#endif 
