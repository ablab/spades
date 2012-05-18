//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "logger/logger.hpp"
#include <simple_tools.hpp>
#include <io/easy_reader.hpp>
#include <io/single_read.hpp>

namespace debruijn_graph {

typedef io::IReader<io::SingleRead> ReadStream;
typedef io::IReader<io::PairedRead> PairedReadStream;
typedef io::MultifileReader<io::PairedRead> MultiPairedStream;
typedef io::MultifileReader<io::SingleRead> MultiSingleStream;

auto_ptr<PairedReadStream> paired_easy_reader(bool followed_by_rc,
		size_t insert_size,
		bool change_read_order = false,
		bool revert_second = true,
		bool original = false,
		io::OffsetType offset_type = io::PhredOffset) {
	vector<PairedReadStream*> streams;
	auto& paired_reads = original ? cfg::get().ds.original_paired_reads : cfg::get().ds.paired_reads;
	for (auto it = paired_reads.begin(); it != paired_reads.end(); ++it) {
		vector<string> filenames = *it;
		io::PairedEasyReader* reader;
		if (filenames.size() == 1) {
			reader = new io::PairedEasyReader(input_file(filenames[0]), followed_by_rc, insert_size, change_read_order, revert_second, offset_type);
		} else if (filenames.size() == 2) {
			auto files = make_pair(input_file((*it)[0]), input_file(filenames[1]));
			reader = new io::PairedEasyReader(files, followed_by_rc, insert_size, change_read_order, revert_second, offset_type);
		} else {
			VERIFY_MSG(false, "Can't handle the case with " << filenames.size() << " input files with paired reads");
		}
		streams.push_back(reader);
	}
	return auto_ptr<PairedReadStream>(new MultiPairedStream(streams, true));
}

auto_ptr<ReadStream> single_easy_reader(bool followed_by_rc,
		bool including_paired_reads,
		bool original = false,
		io::OffsetType offset_type = io::PhredOffset) {
	vector<ReadStream*> streams;
	auto& single_reads = original ? cfg::get().ds.original_single_reads : cfg::get().ds.single_reads;
	if (including_paired_reads) {
		auto& paired_reads = original ? cfg::get().ds.original_paired_reads : cfg::get().ds.paired_reads;
		for (auto it = paired_reads.begin(); it != paired_reads.end(); ++it) {
			for (auto it2 = it->begin(); it2 != it->end(); ++it2) {
				streams.push_back(new io::EasyReader(input_file(*it2), followed_by_rc, offset_type));
				DEBUG("Using input file: " << input_file(*it2));
			}
		}
	}
	for (auto it = single_reads.begin(); it != single_reads.end(); ++it) {
		streams.push_back(new io::EasyReader(input_file(*it), followed_by_rc, offset_type));
		DEBUG("Using input file: " << input_file(*it));
	}
	return auto_ptr<ReadStream>(new MultiSingleStream(streams, true));
}

}
