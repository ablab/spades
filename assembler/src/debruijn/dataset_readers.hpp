/*
 * copy_file.hpp
 *
 *  Created on: 8 Sep 2011
 *      Author: valery
 */
#pragma once

#include <vector>

namespace debruijn_graph {

typedef io::IReader<io::SingleRead> ReadStream;
typedef io::IReader<io::PairedRead> PairedReadStream;
typedef io::MultifileReader<io::SingleRead> MultiFileStream;

auto_ptr<io::PairedEasyReader> paired_easy_reader(size_t insert_size) {
	return auto_ptr<io::PairedEasyReader>(
			new io::PairedEasyReader(
					make_pair(
						input_file(cfg::get().ds.paired_reads[0]),
						input_file(cfg::get().ds.paired_reads[1])
					), insert_size));
}

vector<ReadStream*> single_streams() {
	vector<ReadStream*> streams;
	for (auto it = cfg::get().ds.paired_reads.begin(); it != cfg::get().ds.paired_reads.end(); ++it) {
		streams.push_back(new io::EasyReader(input_file(*it)));
	}
	for (auto it = cfg::get().ds.single_reads.begin(); it != cfg::get().ds.single_reads.end(); ++it) {
		streams.push_back(new io::EasyReader(input_file(*it)));
	}
	return streams;
}

}
