//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "cpp_utils.hpp"

namespace omnigraph {

double get_median(const map<int, size_t> &hist) {
	double S = 0;
	for (auto iter = hist.begin(); iter != hist.end(); ++iter) {
		S += iter->second;
	}
	double sum = S;
	for (auto iter = hist.begin(); iter != hist.end(); ++iter) {
		sum -= iter->second;
		if (sum <= S / 2) {
			return iter->first;
		}
	}
	assert(false);
	return -1;
}

double get_mad(const map<int, size_t> &hist, double median) { // median absolute deviation
	map<int, size_t> hist2;
	for (auto iter = hist.begin(); iter != hist.end(); ++iter) {
		double x = fabs(iter->first - median);
		hist2[x] = iter->second;
	}
	return get_median(hist2);
}

void hist_crop(const map<int, size_t> &hist, double low, double high, map<int, size_t> *res) {
	for (auto iter = hist.begin(); iter != hist.end(); ++iter) {
		if (iter->first >= low && iter->first <= high) {
			res->insert(*iter);
		}
	}
}

template<class graph_pack>
void refine_insert_size(io::IReader<io::PairedRead>& stream, graph_pack& gp, size_t edge_length_threshold) {
	enum {
		k = graph_pack::k_value
	};
	INFO("SUBSTAGE == Refining insert size and its distribution");
	map<int, size_t> hist;
	size_t n = 0;
	size_t total = 0;
	double sum = 0;
	double sum2 = 0;
	INFO("Processing paired reads (takes a while)");
	while (!stream.eof()) {
		io::PairedRead r;
		stream >> r;
		Sequence sequence_left = r.first().sequence();
		Sequence sequence_right = r.second().sequence();
		VERBOSE_POWER(++total, " paired reads processed");
		if (sequence_left.size() <= k || sequence_right.size() <= k) {
			continue;
		}
		Seq<k + 1> left = sequence_left.end<k + 1>();
		Seq<k + 1> right = sequence_right.start<k + 1>();
		left = gp.kmer_mapper.Substitute(left);
		right = gp.kmer_mapper.Substitute(right);
		if (!gp.index.contains(left) || !gp.index.contains(right)) {
			continue; // TODO rather use binary search.
		}
		auto pos_left = gp.index.get(left);
		auto pos_right = gp.index.get(right);
		if (pos_left.first != pos_right.first || gp.g.length(pos_left.first) < edge_length_threshold) {
			continue;
		}
		int is = pos_right.second - pos_left.second - k - 1 - r.insert_size() + sequence_left.size() + sequence_right.size();
		hist[is] += 1;
		++n;
	}

		if (n == 0) {
		throw std::runtime_error("Failed to estimate insert size of paired reads, because none of the paired reads aligned to long edges");
	}
	INFO(n << " paired reads (" << (n * 100.0 / total) << "% of all) aligned to long edges");

	// Misha's approach

	size_t often = 0;
	size_t mode = -1;
	for (auto iter = hist.begin(); iter != hist.end(); ++iter) {
		if (iter->second > often) {
			often = iter->second;
			mode = iter->first;
		}
	}

	int low = -mode;
	int high = 3 * mode;

	n = 0;
	sum = 0;
	sum2 = 0;
	for (auto iter = hist.begin(); iter != hist.end(); ++iter) {
		if (iter->first < low || iter->first > high) {
			continue;
		}
		n += iter->second;
		sum += iter->second * 1.0 * iter->first;
		sum2 += iter->second * 1.0 * iter->first * iter->first;
	}
	double mean, delta;
	mean = sum / n;
	delta = sqrt(sum2 / n - mean * mean);

	low = mean - 5 * delta;
	high = mean + 5 * delta;

	n = 0;
	sum = 0;
	sum2 = 0;
	for (auto iter = hist.begin(); iter != hist.end(); ++iter) {
		if (iter->first < low || iter->first > high) {
			//INFO("outsiders: " << iter->first << " " << iter->second);
			continue;
		}
//		INFO("histogram: " << iter->first << " " << iter->second);
		n += iter->second;
		sum += iter->second * 1.0 * iter->first;
		sum2 += iter->second * 1.0 * iter->first * iter->first;
	}
	mean = sum / n;
	delta = sqrt(sum2 / n - mean * mean);

	size_t m = 0;
	//todo optimize
	size_t q[19];
	for (size_t i = 1; i < 20; ++i) {
		q[i-1] = 5 * i;
	}
	for (auto iter = hist.begin(); iter != hist.end(); ++iter) {
		if (iter->first < low || iter->first > high) {
			continue;
		}
		size_t mm = m + iter->second;
		for (size_t i = 0; i < utils::array_size(q); i++) {
			size_t scaled_q_i(q[i] / 100. * n);
			if (m < scaled_q_i && mm >= scaled_q_i) {
				cfg::get_writable().ds.percentiles[q[i]] = iter->first;
//				INFO("Percentile: " << q[i] << " = " << iter->first);
			}
		}
		m = mm;
	}

	cfg::get_writable().ds.IS = mean;
	cfg::get_writable().ds.is_var = delta;

	// Kolya's approach
	// Now we calculate median, MAD and cropped histogram
	{
		double median = 0;
		double mad = 0;
		median = get_median(hist);
		mad = get_mad(hist, median);
		double low = median - 2 * 1.4826 * mad;
		double high = median + 2 * 1.4826 * mad;
		hist_crop(hist, low, high, &cfg::get_writable().ds.hist);
		median = get_median(cfg::get().ds.hist);
		mad = get_mad(cfg::get().ds.hist, median);
		cfg::get_writable().ds.median = median;
		cfg::get_writable().ds.mad = mad;
	}

	INFO("Insert size refined:");
	INFO("IS = " << cfg::get_writable().ds.IS);
	INFO("delta = " << delta);
}

}
