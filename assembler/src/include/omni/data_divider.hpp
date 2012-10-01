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
    typedef pair<size_t, size_t> Interval;

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
        size_t data_size_;
        int min_value_, max_value_;
        size_t threshold_;

        template <class EdgeId>
        void Print(const vector<PairInfo<EdgeId> >& data) const {
            for (size_t i = 0; i < data.size(); ++i) {
                TRACE(data[i].d << " " << data[i].weight);   
            }
        }

        template <class EdgeId>
        bool IsANewCluster(size_t index, const vector<PairInfo<EdgeId> >& data) {
            VERIFY(index < data_size_ - 1);
            return (abs(data[index + 1].d - data[index].d) > threshold_);
        }

    public:

        DataDivider(size_t threshold) : threshold_(threshold) 
        {
        }

        template <class EdgeId>
        vector<Interval> DivideData(const vector<PairInfo<EdgeId> >& data) {

            data_size_ = data.size();
            VERIFY(data_size_ > 0);
            vector<Interval> answer;

            min_value_ = rounded_d(data.front());
            max_value_ = rounded_d(data.back());
            size_t begin = 0;
            for (size_t i = 0; i < data_size_ - 1; i++) {
                if (IsANewCluster(i, data)){ 
                    answer.push_back(make_pair(begin, i + 1));
                    begin = i + 1;
                }
            }
            answer.push_back(make_pair(begin, data_size_));

            return answer;
        }

        template <class EdgeId>
        vector<Interval> DivideAndSmoothData(const vector<PairInfo<EdgeId> >& data, vector<PairInfo<EdgeId> >& new_data, boost::function<double(int)> weight_f) {

            data_size_ = data.size();
            VERIFY(data_size_ > 0);
            vector<Interval> answer;
            
            TRACE("Data");
            Print(data);

            min_value_ = rounded_d(data.front());
            max_value_ = rounded_d(data.back());
            size_t begin = 0;
            for (size_t i = 0; i < data_size_; i++) {
                if (i == data_size_ - 1 || IsANewCluster(i, data)) {
                    int low_val = data[begin].d;
                    int high_val = data[i].d;
                    size_t new_begin = new_data.size();
                    VERIFY(low_val <= high_val);
                    for (int j = low_val; j <= high_val; ++j) {
                        double val = 0.;
                        for (size_t k = begin; k <= i; ++k) {
                            val += data[k].weight * weight_f(j - data[k].d);
                        }
                        new_data.push_back(PairInfo<EdgeId>(data[0].first, data[0].second, j, val, 0.)); 
                    }
                    size_t new_end = new_data.size();
                    answer.push_back(make_pair(new_begin, new_end));

                    begin = i + 1;
                }
            }
            //answer.push_back(make_pair(beginc, new_data.size()));
            TRACE("New_data ");
            Print(new_data);

            return answer;
        }

    private:
        DECL_LOGGER("DataDivider");
    };
}

#endif /* DATA_DIVIDER_HPP_ */
