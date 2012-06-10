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

//typedef double[2] complex;

template <class EdgeId>
class PeakFinder {

private:

	double x1, x2, y1, y2;

    size_t default_range;
    
    size_t delta_;
	
    double percentage_;

	double derivative_threshold_;

    double weight;

	vector<int> x_;
    vector<double> y_;

	vector<int> peaks;

	int data_size, data_length, min, max;


    typedef std::complex<double> complex_t;

    vector<complex_t> in, out;
	//complex_t *in, *out, *outf;

    int Rev(int num, int lg_n) {
        int res = 0;
        for (int i = 0; i < lg_n; ++i)
            if (num & (1 << i))
                res |= 1 << (lg_n - 1 - i);
                return res;
    }
 
    void FFT(vector<complex_t>& vect, bool invert) {
        int n = (int) vect.size();
        int lg_n = 0;
        while ((1 << lg_n) < n)  
            ++lg_n;

        for (int i = 0; i < n; ++i)
            if (i < Rev(i, lg_n))
                swap(vect[i], vect[Rev(i, lg_n)]);

        for (size_t len = 2; len < 1 + n; len <<= 1) {
            double ang = 2 * M_PI / len * (invert ? -1 : 1);
            complex_t wlen(cos(ang), sin(ang));
            for (size_t i = 0; i < n; i += len) {
                complex_t w(1);
                for (size_t j = 0; j < len >> 1; ++j) {
                    complex_t u = vect[i + j];
                    complex_t v = vect[i + j + len >> 1] * w;
                    vect[i + j] = u + v;
                    vect[i + j + len >> 1] = u - v;
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

    void ExtendLinear(vector<complex_t>& in) {
		int ind = 0;
        weight = 0;
		for (int i = 0; i < data_length; ++i) {
			if (ind == data_size - 1)
				in[i] = max;
			else {
                VERIFY(x_[ind + 1] > x_[ind]);
				in[i] = ((i + min - x_[ind]) * y_[ind + 1] + y_[ind] * (x_[ind + 1] - i - min)) / (1. * (x_[ind + 1] - x_[ind]));
			}
            weight += in[i].real();     // filling the array on the fly

			if (ind < data_size && i == x_[ind + 1] - min) 
                ind++;
		}
	}


	void InitBaseline() {
		int Np = (int) (data_length * percentage_);
		if (Np == 0) Np++; // Np <> 0 !!!!

		double mean_beg = 0.;
		double mean_end = 0.;
		for (int i = 0; i < Np; i++) {
			mean_beg += in[i].real();
			mean_end += in[data_length - i - 1].real();
		}
		mean_beg /= 1. * Np;
		mean_end /= 1. * Np;
		//	two points defining the line
		x1 = Np / 2.;
		x2 = data_length - Np / 2.;
		y1 = mean_beg;
		y2 = mean_end;
	}

	void SubtractBaseline() {
		//	subtracting a baseline
		//	it's being constructed like this: the first point is (Np/2; mean of the first percentage of data),
		//	the second point is (data_length - Np/2; mean of the last $percentage of data)
		for (int i = 0; i < data_length; i++) {
			in[i] -= (y1 + (y2 - y1) * (i - x1) / (x2 - x1));
		}
	}

	void AddBaseline() {
		for (int i = 0; i < data_length; i++) {
			in[i] /= data_length;
			in[i] += (y1 + (y2 - y1) * (i - x1) / (x2 - x1));
		}
	}

	void Init() {
		data_size = x_.size();
		min = x_[0];
		max = x_[data_size - 1];
		data_length = max - min + 1;
		//fast Fourier transform
		//out = (complex_t*) malloc(sizeof(complex_t) * data_length);
		//in = (complex_t*) malloc(sizeof(complex_t) * data_length);
		//outf = (complex_t*) malloc(sizeof(complex_t) * data_length);
		//p = fftw_plan_dft_1d(data_length, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	}

	bool isInRange(int peak) {
		return peak <= max && peak >= min;
	}

	double LeftDerivative(int dist) {
		VERIFY(dist > min);
		return in[dist - min].real() - in[dist - min - 1].real();
	}

	double RightDerivative(int dist) {
		VERIFY(dist < max);
		return in[dist - min + 1].real() - in[dist - min].real();
	}

	double MiddleDerivative(int dist) {
		VERIFY(dist > min && dist < max);
		return 0.5f * (in[dist - min + 1].real() - in[dist - min - 1].real());
	}

	double Derivative(int dist) {
		if (dist == max)
			return LeftDerivative(dist);
		else if (dist == min)
			return RightDerivative(dist);
		else 
            return MiddleDerivative(dist);
    }

	bool isLocalMaximum(int peak, int range, int left_bound, int right_bound, int delta) {

		for (int i = left_bound + range; i <= right_bound - range; i++) {
			int index_max = i - range;
			for (int j = i - range; j <= i + range; j++)
				if (math::ls(in[index_max - min].real(), in[j - min].real())) {
					index_max = j;
				}// else if (j < i && in[index_max - min][0] == in[j - min][0] ) index_max = j;
			

            if (!((index_max > i - (range >> 1)) && (index_max < i + (range >> 1)))) continue;
            if  (abs(index_max - peak) <= delta) return true;
		}
		return false;
    }

	bool isLocalMaximum(int peak, size_t range) {
        return isLocalMaximum(peak, range, min, max, delta_);
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

    double getWeight(){
        return weight;
    }

    double getNormalizedWeight(){
        return weight / data_length;
    }

	void PrintStats() {
		for (int i = 0; i < data_length; i++)
			cout << in[i] << endl;

		cout << endl;
	}

	void FFTSmoothing(double cutoff) {
		VERIFY(data_length > 0);
        if (data_length == 1) {
			in[0] = x_[0];
			in[0] = y_[0];
		}
		ExtendLinear(in);
		InitBaseline();
		SubtractBaseline();
        FFTForward(in);
		//fftw_execute(p);
		//p1 = fftw_plan_dft_1d(data_length, out, outf, FFTW_BACKWARD, FFTW_ESTIMATE);

		int Ncrit = (int) (cutoff);

//      cutting off - standard parabolic filter
        for (int i = 0; i < min(data_length, Ncrit); i++) {
			in[i] *= 1. - (i * i * 1.) / (Ncrit * Ncrit);
		}
		for (int i = Ncrit; i < data_length; i++) {
			in[i] = 0.;
		}

        FFTBackward(in);


		//fftw_execute(p1);
		AddBaseline();
	}

	bool isPeak(int dist, int range) {
        return isLocalMaximum(dist, range = default_range);
	}


//  not tested at all
    vector<pair<int, double> > ListPeaks(int delta = 5) {
        vector<pair<int, double> > peaks;
        //another data_length
        int data_length = max - min + 1;
        bool* was;
        srand(time(NULL));    
        std::fill(was, was + data_length, false);
        for (int l = 0; l < data_length; l++) {
            int v = std::rand() % data_length;
            if (was[v]) continue;
            was[v] = true;
            int index = v + min;
            while (index < max && index > min){
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

            if (index >= max - delta || index <= min + delta) 
                continue;

            if ((right_derivative < -derivative_threshold_) && (left_derivative > derivative_threshold_))
                    if (isLocalMaximum(index, delta, index - delta - 1, index + delta - 1, delta >> 1)) {
                        double weight = 0;
                        for (int i = max(min, index - (delta << 1)); i < min(max, index + (delta << 1)); ++i) {
                            double right_derivative = RightDerivative(i);
                            weight += right_derivative*right_derivative;
                        }
                        peaks.push_back(make_pair(index, weight));
                    }
        }
		return peaks;
	}

	vector<complex_t> getIn() const {
		return in;
	}

	vector<complex_t> getOut() const {
		return in;
	}

	vector<int> getPeaks() const {
		return peaks;
	}

};

}

#endif /* PEAKFINDER_HPP_ */
