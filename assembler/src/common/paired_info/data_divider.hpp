//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
* data_divider.hpp
*
*  Created on: Aug 16, 2011
*      Author: alexeyka
*/


#ifndef DATA_DIVIDER_HPP_
#define DATA_DIVIDER_HPP_

#include <iostream>
#include <math.h>
#include "utils/verify.hpp"
#include <vector>
#include <utility>
#include <cstdlib>
#include <cstdio>
#include "index_point.hpp"

namespace omnigraph {

namespace de {

template<class EdgeId>
class DataDivider {
    typedef pair<size_t, size_t> Interval;
    typedef vector<PairInfo<EdgeId> > PairInfos;
    typedef pair<EdgeId, EdgeId> EdgePair;
    typedef vector<Point> PointArray;
    typedef std::function<double(int)> WeightFunction;

    //    double LeftDerivative(int index, vector<int> x, vector<int> y) {
    //        return outf[dist - min_value_ + 1][0] - outf[dist - min][0];
    //    }
    //
    //    double RightDerivative(index, std::vector<int> x, std::vector<int> y) {
    //        return outf[dist - min_value_][0] - outf[dist - min - 1][0];
    //    }
    //
    //    double MiddleDerivative(int index, std::vector<int> x, std::vector<int> y) {
    //        return 0.5f * (outf[dist - min_value_ + 1][0] - outf[dist - min - 1][0]);
    //    }

public:
    DataDivider(size_t threshold, const PointArray &points) :
            threshold_(threshold), points_(points) {
    }

    vector<Interval> DivideData() {
        VERIFY(points_.size() > 0);
        vector<Interval> answer;
        min_value_ = rounded_d(points_.front());
        max_value_ = rounded_d(points_.back());
        size_t begin = 0;
        for (size_t i = 0; i < points_.size() - 1; ++i) {
            if (IsANewCluster(i, points_)) {
                answer.push_back(make_pair(begin, i + 1));
                begin = i + 1;
            }
        }
        answer.push_back(make_pair(begin, points_.size()));

        return answer;
    }

    vector<Interval> DivideAndSmoothData(const EdgePair &ep,
                                         PairInfos &new_data,
                                         WeightFunction weight_f) {
        VERIFY(points_.size() > 0);
        vector<Interval> answer;

        TRACE("Data");
        //Print();
        const Point &point = points_.front();
        min_value_ = rounded_d(point);
        max_value_ = rounded_d(points_.back());
        size_t begin = 0;
        for (size_t i = 0; i < points_.size(); ++i) {
            if (i == points_.size() - 1 || IsANewCluster(i)) {
                int low_val = rounded_d(points_[begin]);
                int high_val = rounded_d(points_[i]);
                size_t new_begin = new_data.size();
                VERIFY(low_val <= high_val);
                for (int j = low_val; j <= high_val; ++j) {
                    double val = 0.;
                    for (size_t k = begin; k <= i; ++k) {
                        val += points_[k].weight * weight_f(j - rounded_d(points_[k]));
                    }
                    new_data.push_back(PairInfo<EdgeId>(ep.first, ep.second, j, val, 0.));
                }
                size_t new_end = new_data.size();
                answer.push_back(make_pair(new_begin, new_end));

                begin = i + 1;
            }
        }
        //answer.push_back(make_pair(beginc, new_data.size()));
        TRACE("New_data ");
        Print();

        return answer;
    }

private:
    int min_value_;
    int max_value_;
    size_t threshold_;
    PointArray points_;

    void Print() const {
        for (size_t i = 0; i < points_.size(); ++i) {
            TRACE(points_[i].d << " " << points_[i].weight);
        }
    }

    bool IsANewCluster(size_t index) {
        VERIFY(index < points_.size() - 1);
        return (math::gr(abs(points_[index + 1].d - points_[index].d), (DEDistance) threshold_));
    }

    DECL_LOGGER("DataDivider");
};

}


}

#endif /* DATA_DIVIDER_HPP_ */
