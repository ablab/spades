#pragma once

#include "cpp_utils.hpp"

namespace omnigraph {

template<class graph_pack>
void refine_insert_size(io::IReader<io::PairedRead>& stream, graph_pack& gp, size_t edge_length_threshold) {
	enum {
		k = graph_pack::k_value
	};
	INFO("SUBSTAGE == Refining insert size and its distribution");
	map<int, size_t> hist;
	size_t n = 0;
	double sum = 0;
	double sum2 = 0;
	size_t succesfully_processed = 0;
	INFO("Processing paired reads (takes a while)");
	while (!stream.eof()) {
		io::PairedRead r;
		stream >> r;
		Sequence sequence_left = r.first().sequence();
		Sequence sequence_right = r.second().sequence();
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
		succesfully_processed++;
		int is = pos_right.second - pos_left.second - k - 1 - r.insert_size() + sequence_left.size() + sequence_right.size();
		hist[is] += 1;
		VERBOSE_POWER(++n, " paired reads processed");
		sum += is;
		sum2 += is * 1.0 * is;
	}

	const size_t magic_number = 100;
	if (succesfully_processed < magic_number) {
		throw std::runtime_error("Sorry! Failed to estimate paired parameters.");
	}

	double mean, delta;
	mean = sum / n;
	delta = sqrt(sum2 / n - mean * mean);
	size_t often = 0;
	size_t median = -1;
	for (auto iter = hist.begin(); iter != hist.end(); ++iter) {
		if (iter->second > often) {
			often = iter->second;
			median = iter->first;
		}
	}

	int low = -median;
	int high = 3 * median;

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

	INFO("Insert size refined:");
	INFO("IS = " << mean);
	INFO("delta = " << delta);
}

}
