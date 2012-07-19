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
namespace omnigraph {


class DataDivider {
    typedef pair<int, int> Interval;

    //	double LeftDerivative(int index, vector<int> x, vector<int> y) {
    //		return outf[dist - min_value_ + 1][0] - outf[dist - min][0];
    //	}
    //
    //	double RightDerivative(index, std::vector<int> x, std::vector<int> y) {
    //		return outf[dist - min_value_][0] - outf[dist - min - 1][0];
    //	}
    //
    //	double MiddleDerivative(int index, std::vector<int> x, std::vector<int> y) {
    //		return 0.5f * (outf[dist - min_value_ + 1][0] - outf[dist - min - 1][0]);
    //	}

    private:
        int data_size_, data_length_, min_value_, max_value_;
        size_t threshold_;

        bool isCluster(int index, const vector<int> & x, const vector<int> & y) {
            VERIFY(index < data_size_ - 1);
            return (size_t(abs(x[index + 1] - x[index])) > threshold_);
        }

        template <class EdgeId>
        bool isCluster(int index, const vector<PairInfo<EdgeId> >& data) {
            VERIFY(index < data_size_ - 1);
            return (abs(data[index + 1].d - data[index].d) > threshold_);
        }

    public:

        DataDivider(size_t threshold) : threshold_(threshold) 
        {
        }

        template <class EdgeId>
        vector<Interval> DivideData(vector<PairInfo<EdgeId> > data){

            data_size_ = data.size();
            vector<Interval> answer;

            min_value_ = rounded_d(data.front());
            max_value_ = rounded_d(data.back());
            int begin = 0;
            for (int i = 0; i < data_size_ - 1; i++) {
                if (isCluster(i, data)){ 
                    answer.push_back(make_pair(begin, i + 1));
                    begin = i + 1;
                }
            }
            answer.push_back(make_pair(begin, data_size_));

            return answer;
        }

    };
}

#endif /* DATA_DIVIDER_HPP_ */
