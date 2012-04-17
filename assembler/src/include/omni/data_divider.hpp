//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * data_divider.hpp
 *
 *  Created on: Aug 16, 2011
 *      Author: alexeyka
 */

//#include <fftw3.h>
#include <iostream>
#include <math.h>
#include "verify.hpp"
#include <vector>
#include <utility>
#include <cstdlib>
#include <cstdio>
#include "paired_info.hpp"
#include "omni_utils.hpp"

#ifndef DATA_DIVIDER_HPP_
#define DATA_DIVIDER_HPP_
//todo Talk with Alex about this file!!!
namespace omnigraph {

typedef std::pair<int, int> interval;

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

int data_size, data_length, min, max;

size_t Threshold;

bool isCluster(int index, std::vector<int> & x, std::vector<int> & y) {
	VERIFY(index < data_size - 1);
	return (size_t(abs(x[index + 1] - x[index])) > Threshold);
}

template <class EdgeId>
bool isCluster(int index, std::vector<PairInfo<EdgeId> > data) {
	VERIFY(index < data_size - 1);
	return (abs(data[index + 1].d - data[index].d) > Threshold);
}

//static void debug(mydata vec) {
//
//	for (size_t i = 0; i < vec.size(); i++) {
//		std::cout << vec[i].first << " " << vec[i].second << std::endl;
//	}
//	std::cout << std::endl << std::endl;
//
//}

template <class EdgeId>
std::vector<interval> divideData(std::vector<PairInfo<EdgeId> > data){

    data_size = data.size();
    std::vector<interval> answer;

    min = rounded_d(data.front());
    max = rounded_d(data.back());
//	std::cout << "Data size is " << data_size << std::endl;

//	std::cout << "Data length is " << data_length << std::endl;

//	int data_positive_size = data_size >> 1;
//	int data_positive_length = max - data[data_positive_size].d + 1;

//	Thr = (4 * data_positive_length / data_positive_size);

//	std::cout << "Threshold is " << Thr << std::endl << std::endl;
    

    int begin = 0;
    for (int i = 0; i < data_size - 1; i++) {
		if (isCluster(i, data)){ 
            answer.push_back(std::make_pair(begin, i + 1));
            begin = i + 1;
        }
    }
    answer.push_back(std::make_pair(begin, data_size));

	return answer;
}

}

#endif /* DATA_DIVIDER_HPP_ */
