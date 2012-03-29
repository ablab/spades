#pragma once

#include <vector>
#include <logging.hpp>
#include <simple_tools.hpp>

namespace debruijn_graph {

typedef io::IReader<io::SingleRead> ReadStream;
typedef io::IReader<io::PairedRead> PairedReadStream;
typedef io::MultifileReader<io::SingleRead> MultiFileStream;

auto_ptr<io::PairedEasyReader> paired_easy_reader(size_t insert_size) {
	size_t files_count = cfg::get().ds.paired_reads.size();\
	if (files_count == 1) {
		io::PairedEasyReader* reader = new io::PairedEasyReader(input_file(cfg::get().ds.paired_reads[0]), insert_size);
		return auto_ptr<io::PairedEasyReader>(reader);
	}
	if (files_count == 2) {
		auto files = make_pair(input_file(cfg::get().ds.paired_reads[0]), input_file(cfg::get().ds.paired_reads[1]));
		io::PairedEasyReader* reader = new io::PairedEasyReader(files, insert_size);
		return auto_ptr<io::PairedEasyReader>(reader);
	}
	VERIFY_MSG(false, "Can't handle the case with " << ToString(cfg::get().ds.paired_reads.size()) << " input files with paired reads");
}

vector<ReadStream*> single_streams() {
	vector<ReadStream*> streams;
	for (auto it = cfg::get().ds.paired_reads.begin(); it != cfg::get().ds.paired_reads.end(); ++it) {
		streams.push_back(new io::EasyReader(input_file(*it)));
		DEBUG("Using input file: " << input_file(*it));
	}
	for (auto it = cfg::get().ds.single_reads.begin(); it != cfg::get().ds.single_reads.end(); ++it) {
		streams.push_back(new io::EasyReader(input_file(*it)));
		DEBUG("Using input file: " << input_file(*it));
	}
	return streams;
}

}
