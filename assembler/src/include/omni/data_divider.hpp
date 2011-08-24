/*
 * data_divider.hpp
 *
 *  Created on: Aug 16, 2011
 *      Author: alexeyka
 */

#include <fftw3.h>
#include <iostream>
#include <math.h>
#include <assert.h>
#include <vector>
#include <utility>
#include <cstdlib>
#include <cstdio>
#include "paired_info.hpp"
#include "omni_utils.hpp"

#ifndef DATA_DIVIDER_HPP_
#define DATA_DIVIDER_HPP_

namespace omnigraph {

typedef std::pair<int, int> coord;
typedef std::vector<coord> mydata;

//	double LeftDerivative(int index, std::vector<int> x, std::vector<int> y) {
//		return outf[dist - min + 1][0] - outf[dist - min][0];
//	}
//
//	double RightDerivative(index, std::vector<int> x, std::vector<int> y) {
//		return outf[dist - min][0] - outf[dist - min - 1][0];
//	}
//
//	double MiddleDerivative(int index, std::vector<int> x, std::vector<int> y) {
//		return 0.5f * (outf[dist - min + 1][0] - outf[dist - min - 1][0]);
//	}

static int data_size, data_length, min, max;

static double Thr;

static bool isCluster(int index, std::vector<int> & x, std::vector<int> & y) {
	assert(index < data_size - 1);
	return (abs(x[index + 1] - x[index]) > Thr);
}

static bool isCluster(int index, std::vector<PairInfo<EdgeId> > data) {
	assert(index < data_size - 1);
	return (abs(data[index + 1].d - data[index].d) > Thr);
}

//static void debug(mydata vec) {
//
//	for (size_t i = 0; i < vec.size(); i++) {
//		std::cout << vec[i].first << " " << vec[i].second << std::endl;
//	}
//	std::cout << std::endl << std::endl;
//
//}

static std::vector<int> divideData(std::vector<int> x, std::vector<int> y) {
	data_size = x.size();

	std::vector<int> answer;

	answer.push_back(0);
	min = x.front();
	max = x.back();
	data_length = max - min + 1;

	std::cout << "Data size is " << data_size << std::endl;

	std::cout << "Data length is " << data_length << std::endl;

	int data_positive_size = data_size >> 1;
	int data_positive_length = max - x[data_positive_size] + 1;

	Thr = (4 * data_positive_length / data_positive_size);

	std::cout << "Threshold is " << Thr << std::endl << std::endl;

	for (int i = 0; i < data_size - 1; i++) {
		if (isCluster(i, x, y)) answer.push_back(i + 1);
	}

	answer.push_back(data_size);
	return answer;
}

static std::vector<int> divideData(std::vector<PairInfo<EdgeId> > data) {
	data_size = data.size();

	std::vector<int> answer;

	answer.push_back(0);
	min = data.front().d;
	max = data.back().d;
	data_length = max - min + 1;

//	std::cout << "Data size is " << data_size << std::endl;

//	std::cout << "Data length is " << data_length << std::endl;

	int data_positive_size = data_size >> 1;
	int data_positive_length = max - data[data_positive_size].d + 1;

	Thr = (4 * data_positive_length / data_positive_size);

//	std::cout << "Threshold is " << Thr << std::endl << std::endl;

	for (int i = 0; i < data_size - 1; i++) {
		if (isCluster(i, data)) answer.push_back(i + 1);
	}

	answer.push_back(data_size);
	return answer;
}

}

#endif /* DATA_DIVIDER_HPP_ */
