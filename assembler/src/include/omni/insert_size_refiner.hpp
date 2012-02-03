#pragma once

#include "cpp_utils.hpp"

namespace debruijn_graph {

void refine_insert_size(pair<string, string> read_filenames, conj_graph_pack& gp) {
	INFO("SUBSTAGE == Refining insert size and its distribution");
	map<int, size_t> hist;
	size_t n = 0;
	double sum = 0;
	double sum2 = 0;
	io::PairedEasyReader stream(read_filenames,	0);
	while (!stream.eof()) {
		io::PairedRead r;
		stream >> r;
		Sequence sequence_left = r.first().sequence();
		Sequence sequence_right = r.second().sequence();
		if (sequence_left.size() <= K || sequence_right.size() <= K) {
			continue;
		}
		Seq<K + 1> left = sequence_left.end<K + 1>();
		Seq<K + 1> right = sequence_right.start<K + 1>();
		left = gp.kmer_mapper.Substitute(left);
		right = gp.kmer_mapper.Substitute(right);
		if (!gp.index.contains(left) || !gp.index.contains(right)) {
			continue; // TODO rather use binary search.
		}
		auto pos_left = gp.index.get(left);
		auto pos_right = gp.index.get(right);
		if (pos_left.first != pos_right.first) {
			continue;
		}
		int is = pos_right.second - pos_left.second - K - 1 - r.insert_size() + sequence_left.size() + sequence_right.size();
		// DEBUG("insert size refinement evidence: " << is);
		hist[is] += 1;
		n++;
		sum += is;
		sum2 += is * 1.0 * is;
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
		INFO("histogram: " << iter->first << " " << iter->second);
		n += iter->second;
		sum += iter->second * 1.0 * iter->first;
		sum2 += iter->second * 1.0 * iter->first * iter->first;
	}
	mean = sum / n;
	delta = sqrt(sum2 / n - mean * mean);

	size_t m = 0;
	size_t const q[] = {(size_t) (0.05 * n), (size_t) (0.15 * n), (size_t) (0.85 * n), (size_t) (0.95 * n)};
	for (auto iter = hist.begin(); iter != hist.end(); ++iter) {
		if (iter->first < low || iter->first > high) {
			continue;
		}
		size_t mm = m + iter->second;
		for (size_t i = 0; i < utils::array_size(q); i++) {
			if (m < q[i] && mm >= q[i]) {
				INFO("Percentile: " << iter->first);
			}
		}
		m = mm;
	}

	cfg::get_writeable().ds.IS 	= mean;
	INFO("Insert size refined:");
	INFO("IS = " << mean);
	INFO("delta = " << delta);
}

}
