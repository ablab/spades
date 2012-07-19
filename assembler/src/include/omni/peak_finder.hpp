//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * PeakFinder.hpp
 *
 *  Created on: Aug 15, 2011
 *      Author: alexeyka
 */

#ifndef PEAKFINDER_HPP_
#define PEAKFINDER_HPP_

#include "verify.hpp"
#include "data_divider.hpp"
#include "paired_info.hpp"
#include "omni_utils.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <complex>
#include <cmath>


namespace  omnigraph{

template <class EdgeId> class PeakFinder {

private:

	double x1, x2, y1, y2;

    const size_t default_range;
    
    size_t delta_;
	
    double percentage_;

	double derivative_threshold_;

    double weight_;

	vector<int> x_;
    vector<double> y_;

	vector<int> peaks_;

	int data_size_, data_len_, x_left_, x_right_;


    typedef std::complex<double> complex_t;

    vector<complex_t> hist_, out;

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

        for (size_t i = 0; i < n; ++i)
            if (i < Rev(i, lg_n))
                swap(vect[i], vect[Rev(i, lg_n)]);

        for (size_t len = 2; len < 1 + n; len <<= 1) {
            double ang = 2 * M_PI / len * (invert ? -1 : 1);
            complex_t wlen(cos(ang), sin(ang));
            for (size_t i = 0; i < n; i += len) {
                complex_t w(1);
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
                vect[i] /= n;
    }


    void FFTForward(vector<complex_t>& vect) {
       FFT(vect, false);
    }

    void FFTBackward(vector<complex_t>& vect) {
       FFT(vect, true);
    }

    void ExtendLinear(vector<complex_t>& hist) {
		int ind = 0;
        weight_ = 0;
		for (int i = 0; i < data_len_; ++i) {
			if (ind == data_size_ - 1)
				hist.push_back(x_right_);
			else {
                VERIFY(x_[ind + 1] > x_[ind]);
				hist.push_back(((i + x_left_ - x_[ind]) * y_[ind + 1] + y_[ind] * (x_[ind + 1] - i - x_left_)) / (1. * (x_[ind + 1] - x_[ind])));
			}
            weight_ += hist[i].real();     // filling the array on the fly

			if (ind < data_size_ && i == x_[ind + 1] - x_left_) 
                ind++;
		}
	}


	void InitBaseline() {
		int Np = (int) (data_len_ * percentage_);
		if (Np == 0) Np++; // Np <> 0 !!!!

		double mean_beg = 0.;
		double mean_end = 0.;
		for (int i = 0; i < Np; i++) {
			mean_beg += hist_[i].real();
			mean_end += hist_[data_len_ - i - 1].real();
		}
		mean_beg /= 1. * Np;
		mean_end /= 1. * Np;
		//	two points defining the line
		x1 = Np / 2.;
		x2 = data_len_ - Np / 2.;
		y1 = mean_beg;
		y2 = mean_end;
	}

	void SubtractBaseline() {
		//	subtracting a baseline
		//	it's being constructed like this: the first point is (Np/2; mean of the first percentage of data),
		//	the second point is (data_len_ - Np/2; mean of the last $percentage of data)
		for (int i = 0; i < data_len_; i++) {
			hist_[i] -= (y1 + (y2 - y1) * (i - x1) / (x2 - x1));
		}
	}

	void AddBaseline() {
		for (int i = 0; i < data_len_; i++) {
			hist_[i] /= data_len_;
			hist_[i] += (y1 + (y2 - y1) * (i - x1) / (x2 - x1));
		}
	}

	void Init() {
		data_size_ = x_.size();
		x_left_ = x_[0];
		x_right_ = x_[data_size_ - 1];
		data_len_ = x_right_ - x_left_ + 1;
		//fast Fourier transform
		//out = (complex_t*) malloc(sizeof(complex_t) * data_len_);
		//hist_ = (complex_t*) malloc(sizeof(complex_t) * data_len_);
		//outf = (complex_t*) malloc(sizeof(complex_t) * data_len_);
		//p = fftw_plan_dft_1d(data_len_, hist_, out, FFTW_FORWARD, FFTW_ESTIMATE);
	}

	bool isInRange(int peak) {
		return peak <= x_right_ && peak >= x_left_;
	}

	double LeftDerivative(int dist) {
		VERIFY(dist > x_left_);
		return hist_[dist - x_left_].real() - hist_[dist - x_left_ - 1].real();
	}

	double RightDerivative(int dist) {
		VERIFY(dist < x_right_);
		return hist_[dist - x_left_ + 1].real() - hist_[dist - x_left_].real();
	}

	double MiddleDerivative(int dist) {
		VERIFY(dist > x_left_ && dist < x_right_);
		return .5 * (hist_[dist - x_left_ + 1].real() - hist_[dist - x_left_ - 1].real());
	}

	double Derivative(int dist) {
		if (dist == x_right_)
			return LeftDerivative(dist);
		else if (dist == x_left_)
			return RightDerivative(dist);
		else 
            return MiddleDerivative(dist);
    }

	bool isLocalMaximum(int peak, int range, int left_bound, int right_bound, int delta) {

		for (int i = left_bound + range; i <= right_bound - range; i++) {
			int index_max = i - range;
			for (int j = i - range; j <= i + range; j++)
				if (math::ls(hist_[index_max - x_left_].real(), hist_[j - x_left_].real())) {
					index_max = j;
				}// else if (j < i && hist_[index_max - x_left_][0] == hist_[j - x_left][0] ) index_max = j;
			

            if (!((index_max > i - (range >> 1)) && (index_max < i + (range >> 1)))) continue;
            if  (abs(index_max - peak) <= delta) return true;
		}
		return false;
    }

	bool isLocalMaximum(int peak, size_t range) {
        return isLocalMaximum(peak, range, x_left_, x_right_, delta_);
    }

public:
	PeakFinder(vector<PairInfo<EdgeId> > data, int begin, int end, size_t range, size_t delta, double percentage, double derivative_threshold) : 
            default_range(range), 
            delta_(delta), 
            percentage_(percentage), 
            derivative_threshold_(derivative_threshold)
    {
        for (int i = begin; i < end; i++) {
            x_.push_back(rounded_d(data[i]));
            y_.push_back(data[i].weight);
        }
        Init();
    }

    double weight(){
        return weight_;
    }

    double getNormalizedWeight(){
        return weight_ / data_len_;
    }

	void PrintStats() {
		for (int i = 0; i < data_len_; i++)
			cout << hist_[i] << endl;

		cout << endl;
	}

	void FFTSmoothing(double cutoff) {
		VERIFY(data_len_ > 0);
        if (data_len_ == 1) {
			hist_[0] = x_[0];
			hist_[0] = y_[0];
		}
		ExtendLinear(hist_);
		InitBaseline();
		SubtractBaseline();
        FFTForward(hist_);
		int Ncrit = (int) (cutoff);

//      cutting off - standard parabolic filter
        for (int i = 0; i < data_len_ && i < Ncrit; ++i) {
			hist_[i] *= 1. - (i * i * 1.) / (Ncrit * Ncrit);
		}
		for (int i = Ncrit; i < data_len_; i++) {
			hist_[i] = 0.;
		}
        FFTBackward(hist_);
		AddBaseline();
	}

	bool isPeak(int dist, int range) {
        return isLocalMaximum(dist, range);
	}

	bool isPeak(int dist) {
        return isLocalMaximum(dist, default_range);
	}

//  not tested at all
    vector<pair<int, double> > ListPeaks(int delta = 5) {
        vector<pair<int, double> > peaks_;
        //another data_len_
        int data_len_ = x_right_ - x_left_ + 1;
        bool* was;
        srand(time(NULL));    
        std::fill(was, was + data_len_, false);
        for (int l = 0; l < data_len_; l++) {
            int v = std::rand() % data_len_;
            if (was[v]) continue;
            was[v] = true;
            int index = v + x_left_;
            while (index < x_right_ && index > x_left_){
                // if @index is local maximum, then leave it
                double right_derivative = RightDerivative(index);
                double left_derivative = LeftDerivative(index);

                if (right_derivative > 0 && right_derivative >= left_derivative){
                    index++;
                }else if (left_derivative > 0){
                    index--;
                }
            }

            double right_derivative = RightDerivative(index);
            double left_derivative = LeftDerivative(index);

            if (index >= x_right_ - delta || index <= x_left_ + delta) 
                continue;

            if ((right_derivative < -derivative_threshold_) && (left_derivative > derivative_threshold_))
                    if (isLocalMaximum(index, delta, index - delta - 1, index + delta - 1, delta >> 1)) {
                        double weight_ = 0;
                        for (int i = max(x_left_, index - (delta << 1)); i < max(x_right_, index + (delta << 1)); ++i) {
                            double right_derivative = RightDerivative(i);
                            weight_ += right_derivative * right_derivative;
                        }
                        peaks_.push_back(make_pair(index, weight_));
                    }
        }
		return peaks_;
	}

	vector<complex_t> getIn() const {
		return hist_;
	}

	vector<complex_t> getOut() const {
		return hist_;
	}

	vector<int> getPeaks() const {
		return peaks_;
	}

};

}

#endif /* PEAKFINDER_HPP_ */
