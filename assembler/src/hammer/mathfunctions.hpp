//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * mathfunctions.hpp
 *
 *  Created on: 11.05.2011
 *      Author: snikolenko
 */

#ifndef MATHFUNCTIONS_H
#define MATHFUNCTIONS_H

#include "kmer_data.hpp"

#include <cmath>

inline long double logSimplexVolume(int n) {
	return (-lgamma(n+1));
}

inline long double Factorial(int n) {
  if (n == 0) {
    return 1;
  }
  static unordered_map<int, long double> ans;
  if (ans.count(n) == 0) {
    ans[n] = Factorial(n - 1) * n;
  }
  return ans[n];
}

inline long double CNK(int n, int k) {
  return Factorial(n) / (Factorial(k) * Factorial(n - k));
}

inline long double Bernoulli(int k, int n, long double p) {
  return pow(p, k) * pow(1 - p, n - k) * CNK(n, k);
}


/**
  * @return logarithm of {n choose k}
  */
inline double logChooseNK(int n, int k) {
	double res = 0;
	for (int i=1; i <= k; ++i) {
		res += log(n-k+i) - log(i);
	}
	return res;
}

/**
  * @return Beta(x,y)
  */
inline double lBeta(int x, int y) {
	return (lgamma(x) + lgamma(y) - lgamma(x+y));
}

/**
  * @return log({a_1+...+a_n \choose a_1, ..., a_n})
  */
inline double lMultinomial(const vector<unsigned> & x, const KMerData & data_) {
	double res = 0.0, sum = 0.0;
	for (size_t i=0; i<x.size(); ++i) {
		res += lgamma(data_[x[i]].count+1);
		sum += data_[x[i]].count;
	}
	return (lgamma(sum+1) - res);
}

/**
  * @return log({a_1+...+a_n \choose a_1, ..., a_n})
  */
inline double lMultinomial(const vector<KMerStat> & x) {
	double res = 0.0, sum = 0.0;
	for (size_t i=0; i<x.size(); ++i) {
		res += lgamma(x[i].count+1);
		sum += x[i].count;
	}
	return (lgamma(sum+1) - res);
}

/**
  * @return log({a_1+...+a_n \choose a_1, ..., a_n})
  */
inline double lMultinomial(const vector<StringCount> & x) {
	double res = 0.0, sum = 0.0;
	for (size_t i=0; i<x.size(); ++i) {
		res += lgamma(x[i].second.first+1);
		sum += x[i].second.first;
	}
	return (lgamma(sum+1) - res);
}

/**
  * @return log({a_1+...+a_n \choose a_1, ..., a_n}) for reads corresponding to the mask
  */
inline double lMultinomialWithMask(const vector<unsigned> & x, const KMerData &data_, const vector<int> & mask, int maskval) {
	assert(x.size() == mask.size());
	double res = 0.0, sum = 0.0;
	for (size_t i=0; i<x.size(); ++i) {
		if (mask[i] != maskval) continue;
		res += lgamma(data_[x[i]].count+1);
		sum += data_[x[i]].count;
	}
	return (lgamma(sum+1) - res);
}

/**
  * @return log(Beta(a_1+1, ..., a_n+1))
  */
inline double lBetaPlusOne(const vector<int> & x) {
	double res = 0.0, sum = 0.0;
	for (size_t i=0; i<x.size(); ++i) {
		res += lgamma(x[i]+1);
		sum += x[i]+1;
	}
	return (res - lgamma(sum));
}

#endif 
