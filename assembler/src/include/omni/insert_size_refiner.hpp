#pragma once

namespace debruijn_graph {

void refine_insert_size(pair<string, string> read_filenames, conj_graph_pack& gp) {
	if (cfg::get().ds.IS && cfg::get().ds.delta) {
		return;
	}
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
		if (!gp.index.containsInIndex(left) || !gp.index.containsInIndex(right)) {
			continue; // TODO rather use binary search.
		}
		auto pos_left = gp.index.get(left);
		auto pos_right = gp.index.get(right);
		if (pos_left.first != pos_right.first) {
			continue;
		}
		int is = pos_right.second - pos_left.second - K - 1 - r.insert_size() + sequence_left.size() + sequence_right.size();
		TRACE("refine insert size evidence: " << is);
		hist[is] += 1;
		n++;
		sum += is;
		sum2 += is * 1.0 * is;
	}
	double mean = sum / n;
	double delta = sqrt(sum2 / n - mean * mean);
	cfg::get_writeable().ds.IS 	= mean;
	cfg::get_writeable().ds.delta = delta;
	cfg::get_writeable().ds.is_refined = true;
	size_t often = 0;
	size_t median = -1;
	for (auto iter = hist.begin(); iter != hist.end(); ++iter) {
		INFO("histogram: " << iter->first << " " << iter->second);
		if (iter->second > often) {
			often = iter->second;
			median = iter->first;
		}
	}
	n = 0;
	sum = 0;
	sum2 = 0;
	for (auto iter = hist.begin(); iter != hist.end(); ++iter) {
		if (iter->first < -mean || iter->first > 3 * mean) {
			continue;
		}
		n += iter->second;
		sum += iter->second * 1.0 * iter->first;
		sum2 += iter->second * 1.0 * iter->first * iter->first;
	}
	mean = sum / n;
	delta = sqrt(sum2 / n - mean * mean);
	for (auto iter = hist.begin(); iter != hist.end(); ++iter) {
		if (iter->first < mean - 5 * delta || iter->first > mean + 5 * delta) {
			continue;
		}
		n += iter->second;
		sum += iter->second * 1.0 * iter->first;
		sum2 += iter->second * 1.0 * iter->first * iter->first;
	}
	INFO("Insert size refined:");
	INFO("IS = " << *cfg::get().ds.IS);
	INFO("delta = " << *cfg::get().ds.delta);
}

}
