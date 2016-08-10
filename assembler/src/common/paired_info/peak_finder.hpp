//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * PeakFinder.hpp
 *
 *  Created on: Aug 15, 2011
 *      Author: alexeyka
 */

#ifndef PEAKFINDER_HPP_
#define PEAKFINDER_HPP_

#include "utils/verify.hpp"
#include "data_divider.hpp"
#include "paired_info.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <complex>
#include <cmath>

namespace  omnigraph{

namespace de {

template <class EdgeId>
class PeakFinder {

  typedef std::complex<double> complex_t;

 public:
    PeakFinder(const vector<PairInfo<EdgeId> >& data,
             size_t begin,
             size_t end,
             size_t /*range*/,
             size_t delta,
             double percentage,
             double der_thr) :
                delta_(delta),
                percentage_(percentage),
                der_thr_(der_thr)
  {
    for (size_t i = begin; i < end; ++i) {
      x_.push_back(rounded_d(data[i]));
      y_.push_back(data[i].weight());
    }
    Init();
  }

  double weight() const {
    return weight_;
  }

  double GetNormalizedWeight() const {
    return weight_;
  }

  void PrintStats(string host) const {
    for (size_t i = 0; i < data_len_; ++i)
      DEBUG(host << (x_left_ + (int) i) << " " << hist_[i]);
  }

  void FFTSmoothing(double cutoff) {
    VERIFY(data_len_ > 0);
    if (data_len_ == 1) {
      hist_[0] = y_[0];
      return;
    }
    InitBaseline();
    SubtractBaseline();
    FFTForward(hist_);
    size_t Ncrit = (size_t) (cutoff);

    //      cutting off - standard parabolic filter
    for (size_t i = 0; i < data_len_ && i < Ncrit; ++i)
      hist_[i] *= 1. - ((double) i * (double) i * 1.) / (double) (Ncrit * Ncrit);

    for (size_t i = Ncrit; i < hist_.size(); ++i)
      hist_[i] = 0.;

    FFTBackward(hist_);
    AddBaseline();
  }

  bool IsPeak(int dist, size_t range) const {
    return IsLocalMaximum(dist, range);
  }

  bool IsPeak(int dist) const {
    return IsLocalMaximum(dist, 10);
  }

  //  looking for one maximum in the picture
  vector<pair<int, double> > ListPeaks(/*int delta = 3*/) const {
    TRACE("Smoothed data");
    //size_t index_max = 0;
    //for (size_t i = 0; i < data_len_; ++i) {
    //TRACE(x_left_ + (int) i << " " << hist_[i]);
    //if (hist_[i].real() > hist_[index_max].real())
    //index_max = i;
    //}
    //vector<pair<int, double> > result;
    //result.push_back(make_pair(x_left_ + index_max, hist_[index_max].real()));
    //return result;
    DEBUG("Listing peaks");
    map<int, double> peaks_;
    //another data_len_
    size_t data_len_ = (size_t) (x_right_ - x_left_);
    vector<bool> was;
    srand((unsigned) time(NULL));
    for (size_t i = 0; i < data_len_; ++i)
      was.push_back(false);

    size_t iteration = 0;
    for (size_t l = 0; l < data_len_; ++l) {
    //for (size_t k = 0; k < 4; ++k) {
      //size_t v = std::rand() % data_len_;
      size_t v = l;
      if (was[v])
        continue;

      was[v] = true;
      int index = (int) v + x_left_;
      while (index < (x_right_ - 1) && index > x_left_ && iteration < 5) {
        // if @index is local maximum, then leave it
        double right_derivative = RightDerivative(index);
        double left_derivative = LeftDerivative(index);

        if (math::gr(right_derivative, 0.) && math::gr(right_derivative, -left_derivative)) {
          index++;
          if ((iteration & 1) == 0)
            ++iteration;
        }
        else if (math::le(left_derivative, 0.)) {
          index--;
          if ((iteration & 1) == 1)
            ++iteration;
        }
        else
          break;
      }

      TRACE("FOUND " << index);

      //double right_derivative = RightDerivative(index);
      //double left_derivative = LeftDerivative(index);

      if (index < 0)
        continue;

      //if (index >= x_right_ - delta || index < x_left_ + delta)
      //continue;

      TRACE("Is in range");

      if (IsPeak(index, 5)) {
        TRACE("Is local maximum " << index);
        double weight_ = 0.;
        int left_bound = (x_left_ > (index - 20) ? x_left_ : (index - 20));
        int right_bound = (x_right_ < (index + 1 + 20) ? x_right_ : (index + 1 + 20));
        for (int i = left_bound; i < right_bound; ++i)
          weight_ += hist_[i - x_left_].real();
        TRACE("WEIGHT counted");
        pair<int, double> tmp_pair = make_pair(index, 100. * weight_);
        if (!peaks_.count(index)) {
          TRACE("Peaks size " << peaks_.size() << ", inserting " << tmp_pair);
          peaks_.insert(tmp_pair);
        } else {
          TRACE("NON UNIQUE");
        }
      }
    }
    TRACE("FINISHED " << peaks_.size());
    vector<pair<int, double> > peaks;
    for (auto iter = peaks_.begin(); iter != peaks_.end(); ++iter) {
      const pair<int, double>& tmp_pair = *iter;
      TRACE("next peak " << tmp_pair);
      peaks.push_back(tmp_pair);
      //for (int i = -10; i <= 10; ++i) {
      //peaks.push_back(make_pair(tmp_pair.first + i, tmp_pair.second / 21.));
      //}
    }
    return peaks;
  }

    vector<complex_t> getIn() const {
        return hist_;
    }

    vector<complex_t> getOut() const {
        return hist_;
    }

private:
  double x1, x2, y1, y2;
  size_t delta_;
  double percentage_;
  double der_thr_;
  double weight_;
  vector<int> x_;
  vector<double> y_;
  size_t data_size_, data_len_;
  int x_left_, x_right_;
  vector<complex_t> hist_;

  size_t Rev(size_t num, size_t lg_n) {
    size_t res = 0;
    for (size_t i = 0; i < lg_n; ++i)
      if (num & (1 << i))
        res |= 1 << (lg_n - 1 - i);
    return res;
  }

  void FFT(vector<complex_t>& vect, bool invert) {
    size_t n = vect.size();
    size_t lg_n = 0;
    while ( (1u << lg_n) < n)
      ++lg_n;

    while (n < (1u << lg_n)) {
      vect.push_back(0.);
      ++n;
    }

    for (size_t i = 0; i < n; ++i)
      if (i < Rev(i, lg_n))
        swap(vect[i], vect[Rev(i, lg_n)]);

    for (size_t len = 2; len < 1 + n; len <<= 1) {
      double ang = 2 * M_PI / (double) len * (invert ? -1 : 1);
      complex_t wlen(cos(ang), sin(ang));
      for (size_t i = 0; i < n; i += len) {
        complex_t w(1.);
        for (size_t j = 0; j < (len >> 1); ++j) {
          complex_t u = vect[i + j];
          complex_t v = vect[i + j + (len >> 1)] * w;
          vect[i + j] = u + v;
          vect[i + j + (len >> 1)] = u - v;
          w *= wlen;
        }
      }
    }

    if (invert)
      for (size_t i = 0; i < n; ++i)
        vect[i] /= (double) n;
  }


  void FFTForward(vector<complex_t>& vect) {
    FFT(vect, false);
  }

  void FFTBackward(vector<complex_t>& vect) {
    FFT(vect, true);
  }

  void ExtendLinear(vector<complex_t>& hist) {
    size_t ind = 0;
    weight_ = 0.;
    for (size_t i = 0; i < data_len_; ++i) {
      if (ind == data_size_ - 1)
        hist.push_back((double) x_right_);
      else {
        VERIFY(x_[ind + 1] > x_[ind]);
        hist.push_back(((double) (i + x_left_ - x_[ind]) *
                        y_[ind + 1] + y_[ind] *
                        (double) (x_[ind + 1] - i - x_left_)) /
                        (double) (1 * (x_[ind + 1] - x_[ind])));
      }
      weight_ += hist[i].real();     // filling the array on the fly

      if (ind < data_size_ && ((int) i == x_[ind + 1] - x_left_))
        ++ind;
    }

  }


  void InitBaseline() {
    size_t Np = (size_t) ((double) data_len_ * percentage_);
    if (Np == 0) Np++; // Np <> 0 !!!!

    double mean_beg = 0.;
    double mean_end = 0.;
    for (size_t i = 0; i < Np; ++i) {
      mean_beg += hist_[i].real();
      mean_end += hist_[data_len_ - i - 1].real();
    }
    mean_beg /= 1. * (double) Np;
    mean_end /= 1. * (double) Np;

    //    two points defining the line
    x1 = (double) Np / 2.;
    x2 = (double) data_len_ - (double) Np / 2.;
    y1 = mean_beg;
    y2 = mean_end;
  }

  void SubtractBaseline() {
    //    subtracting a baseline
    //    it's being constructed like this: the first point is (Np/2; mean of the first percentage of data),
    //    the second point is (data_len_ - Np/2; mean of the last $percentage of data)
    for (size_t i = 0; i < data_len_; ++i) {
      hist_[i] -= (y1 + (y2 - y1) * ((double) i - x1) / (x2 - x1));
    }
  }

  void AddBaseline() {
    for (size_t i = 0; i < data_len_; ++i) {
      hist_[i] += (y1 + (y2 - y1) * ((double) i - x1) / (x2 - x1));
    }
  }

  void Init() {
    data_size_ = x_.size();
    x_left_ = x_[0];
    x_right_ = x_[data_size_ - 1] + 1;
    data_len_ = x_right_ - x_left_;
    ExtendLinear(hist_);
  }

  bool IsInRange(int peak) const {
    return peak < x_right_ && peak >= x_left_;
  }

  double LeftDerivative(int dist) const {
    VERIFY(dist > x_left_);
    return hist_[dist - x_left_].real() - hist_[dist - x_left_ - 1].real();
  }

  double RightDerivative(int dist) const {
    VERIFY(dist < x_right_ - 1);
    return hist_[dist - x_left_ + 1].real() - hist_[dist - x_left_].real();
  }

  double MiddleDerivative(int dist) const {
    VERIFY(dist > x_left_ && dist < x_right_ - 1);
    return .5 * (hist_[dist - x_left_ + 1].real() - hist_[dist - x_left_ - 1].real());
  }

  double Derivative(int dist) const {
    if (dist == x_right_ - 1)
      return LeftDerivative(dist);
    else if (dist == x_left_)
      return RightDerivative(dist);
    else
      return MiddleDerivative(dist);
  }

  bool IsLocalMaximum(int peak, size_t range, int left_bound, int right_bound, size_t delta) const {

    DEBUG("Is local maximum :  peak " << peak << " range " << range
       << " bounds " << left_bound << " " << right_bound << " delta " << delta);
    int index_max = peak;
    TRACE("Looking for the maximum");
    for (int j = left_bound; j < right_bound; ++j)
      if (math::ls(hist_[index_max - x_left_].real(), hist_[j - x_left_].real())) {
        index_max = j;
      }// else if (j < i && hist_[index_max - x_left_][0] == hist_[j - x_left][0] ) index_max = j;
    TRACE("Maximum is " << index_max);

    if  ((size_t)abs(index_max - peak) <= delta)
      return true;

    return false;
  }

  bool IsLocalMaximum(int peak, size_t range) const {
    return IsLocalMaximum(peak, range, x_left_, x_right_, delta_);
  }

  DECL_LOGGER("PeakFinder");
};

}

}

#endif /* PEAKFINDER_HPP_ */
