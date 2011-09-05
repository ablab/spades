/*
 * PeakFinder.hpp
 *
 *  Created on: Aug 15, 2011
 *      Author: alexeyka
 */

#ifndef PEAKFINDER_HPP_
#define PEAKFINDER_HPP_

#include <fftw3.h>
#include <assert.h>
#include "data_divider.hpp"
#include "paired_info.hpp"
#include "omni_utils.hpp"



namespace  omnigraph{


typedef fftw_complex complex;

class PeakFinder {
private:

	double x1, x2, y1, y2;

	static const double Percentage;

	static const double DerivativeThreshold;

	static const double HeightThreshold;

//	static const int range;

	std::vector<int> x_, y_;

	std::vector<int> peaks;

	int data_size, data_length, min, max;

	complex *in, *out, *outf;

	fftw_plan p, p1;

	inline void ExtendLinear(complex* f) {
		int ind = 0;
		for (int i = 0; i < data_length; i++) {
			if (ind == data_size - 1)
				f[i][0] = max;
			else {
				f[i][0] = ((i + min - x_[ind]) * y_[ind + 1] + y_[ind] * (x_[ind + 1] - i - min)) / (1.0f * (x_[ind + 1] - x_[ind]));
			}
			f[i][1] = 0;
//			cout<<f[i][0]<<endl;

			if (ind < data_size && i == x_[ind + 1] - min) ind++;
		}
	}

	inline void ExtendNaive(complex* f) {
		int min = x_.front();
		int data_size = x_.size();
		for (int i = 0; i < data_size; i++) {
			f[-min + x_[i]][0] = y_[i];
			f[-min + x_[i]][1] = 0;
			if (i < data_size - 1) for (int j = x_[i] - min + 1; j < x_[i + 1] - min; j++)
				f[j][0] = f[j][1] = 0;
		}
//		for (int i = 0; i < data_length; i++)
//			std::cout << f[i][0] << std::endl;
	}

	void InitBaseline() {
		int Np = (int) ((((data_length * Percentage))));
		if (Np == 0) Np++; // Np <> 0 !!!!

		double mean_beg = 0.0;
		double mean_end = 0.0;
		for (int i = 0; i < Np; i++) {
			mean_beg += in[i][0];
			mean_end += in[data_length - i - 1][0];
		}
		mean_beg /= 1.0f * Np;
		mean_end /= 1.0f * Np;
		//	two points defining the line
		x1 = Np / 2.0f;
		x2 = data_length - Np / 2.0f;
		y1 = mean_beg;
		y2 = mean_end;
	}

	void SubtractBaseline() {
		//	subtracting a baseline
		//	it's being constructed like this: the first point is (Np/2; mean of the first percentage of data),
		//	the second point is (data_length - Np/2; mean of the last $percentage of data)
		for (int i = 0; i < data_length; i++) {
			in[i][0] -= (y1 + (y2 - y1) * (i - x1) / (x2 - x1));
		}
	}

	void AddBaseline() {
		for (int i = 0; i < data_length; i++) {
			outf[i][0] /= data_length;
			outf[i][0] += (y1 + (y2 - y1) * (i - x1) / (x2 - x1));
		}
	}

	void Init() {
		data_size = x_.size();
		min = x_[0];
		max = x_[data_size - 1];
		data_length = max - min + 1;
//		std::cout << data_size << " " << data_length << std::endl;
		//fast Fourier transform
		out = (fftw_complex*) (((fftw_malloc(sizeof(fftw_complex) * data_length))));
		in = (fftw_complex*) (((fftw_malloc(sizeof(fftw_complex) * data_length))));
		outf = (fftw_complex*) (((fftw_malloc(sizeof(fftw_complex) * data_length))));
		p = fftw_plan_dft_1d(data_length, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	}

	bool isInrange(int peak) {
		return peak <= max && peak >= min;
	}

	double LeftDerivative(int dist) {
		assert(dist>min);
		return outf[dist - min][0] - outf[dist - min - 1][0];
	}

	double RightDerivative(int dist) {
		assert(dist<max);
		return outf[dist - min + 1][0] - outf[dist - min][0];
	}

	double MiddleDerivative(int dist) {
		assert(dist>min && dist<max);
		return 0.5f * (outf[dist - min + 1][0] - outf[dist - min - 1][0]);
	}

	double Derivative(int dist) {
		if (dist == max)
			return LeftDerivative(dist);

		else if (dist == min)
			return RightDerivative(dist);

		else return MiddleDerivative(dist);

	}
	bool isLocalMaximum(int peak, int range) {
//		int left_limit = std::max(peak - range * 2, min);
//		int right_limit = std::min(max, peak + range * 2);

		for (int i = min + range; i <= max - range; i++) {
			int index_max = i - range;
			for (int j = i - range; j <= i + range; j++)
				if (outf[index_max - min][0] < outf[j - min][0]) {
					index_max = j;
				} else if (j < i && outf[index_max - min][0] == outf[j - min][0] ) index_max = j;
			

            if (!((index_max > i - (range >> 1)) && (index_max < i + (range >> 1)))) continue;
//            if (RightDerivative(index_max + 1)>-DerivativeThreshold || LeftDerivative(index_max - 1)<DerivativeThreshold) 
//            continue;
//            std::cout<< RightDerivative(index_max + 1) << " HUISHUIS "<< LeftDerivative(index_max - 1) << std::endl;
            if  (abs(index_max - peak) < data_length>>1) return true;
		}
		return false;
    
    }

public:
	PeakFinder(std::vector<int> x, std::vector<int> y) :
			x_(x), y_(y) {
		Init();
	}

	PeakFinder(std::vector<int> x, std::vector<int> y, int begin, int end) {
//		std::cout << "hahaha " << begin << " " << end << std::endl;
		for (int i = begin; i < end; i++) {
			x_.push_back(x[i]);
			y_.push_back(y[i]);
		}
		Init();
	}

	PeakFinder(std::vector<PairInfo<EdgeId> > data, int begin, int end) {
	//		std::cout << "hahaha " << begin << " " << end << std::endl;
			for (int i = begin; i < end; i++) {
				x_.push_back(rounded_d(data[i]));
				y_.push_back(data[i].weight);
			}
			Init();
		}
	~PeakFinder() {
		fftw_destroy_plan(p);
		fftw_destroy_plan(p1);
		fftw_free(in);
		fftw_free(out);
		fftw_free(outf);
	}

	void printstat() {
		for (int i = 0; i < data_length; i++)
			std::cout << in[i][0] << " " << in[i][1] << "*I" << std::endl;

		std::cout << std::endl;
	}
	void FFTSmoothing(double cutoff) {
		assert(data_length>0);
        if (data_length == 1){
			in[0][0] = x_[0];
			outf[0][0] = y_[0];
		}
		// linear extension right HERE, because it'll be faster
		ExtendLinear(in);
		InitBaseline();
		SubtractBaseline();
		fftw_execute(p);
		//	cutting off
		p1 = fftw_plan_dft_1d(data_length, out, outf, FFTW_BACKWARD, FFTW_ESTIMATE);

		int Ncrit = (int) (cutoff);
//		cout<< "NCRITICAL "<<Ncrit << std::endl;
		for (int i = 0; i < std::min(data_length, Ncrit); i++) {
			out[i][0] *= 1 - (i * i * 1.0f) / (Ncrit * Ncrit);
			out[i][1] *= 1 - (i * i * 1.0f) / (Ncrit * Ncrit);
		}
		for (int i = Ncrit; i < data_length; i++) {
			out[i][0] = out[i][1] = 0;
		}


		fftw_execute(p1);
		AddBaseline();
		//		for (int i = 0; i < data_length; i+Ф+)
		//			std::cout << in[i][0] << " " << in[i][1] << "*I" << "           " << out[i][0] << " " << out[i][1] << "*I" << "           "
		//					<< out1[i][0] / data_length << " " << out1[i][1] / data_length << "*I" << std::endl;
	}

	bool isPeak(int dist) {
		if (!isInrange(dist)) return false;

        int range = data_length >> 2;

        if (isLocalMaximum(dist, range)) return true;

		return false;
	}
	std::vector<int> ListPeaks() {
		for (int i = min; i <= max; i++) {
			if (isPeak(i)) {
				peaks.push_back(i);
//				cout << "PEEEEEEEEAK : " << i << endl;
			}
		}

		return peaks;
	}

	complex *getIn() const {
		return in;
	}

	void setIn(complex *in) {
		this->in = in;
	}

	complex *getOut() const {
		return outf;
	}

	std::vector<int> getPeaks() const {
		return peaks;
	}

};

const double PeakFinder::Percentage = 0.01f;

const double PeakFinder::DerivativeThreshold = 0.1f;

const double PeakFinder::HeightThreshold = 1;
}

#endif /* PEAKFINDER_HPP_ */
