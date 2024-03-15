//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "library_data.hpp"

#include "llvm/Support/YAMLTraits.h"

using namespace debruijn_graph::config;

namespace llvm { namespace yaml {

void MappingTraits<LibraryData::BinaryReadsInfo>::mapping(IO &io, LibraryData::BinaryReadsInfo &info) {
    io.mapRequired("binary converted", info.binary_converted);
    io.mapRequired("bin reads info file", info.bin_reads_info_file);
    io.mapRequired("paired read prefix", info.paired_read_prefix);
    io.mapRequired("merged read prefix", info.merged_read_prefix);
    io.mapRequired("single read prefix", info.single_read_prefix);
    io.mapRequired("chunk num", info.chunk_num);
}

void MappingTraits<LibraryData>::mapping(IO &io, debruijn_graph::config::LibraryData &data) {
    io.mapRequired("unmerged read length", data.unmerged_read_length);
    io.mapRequired("merged read length", data.merged_read_length);
    io.mapRequired("insert size mean", data.mean_insert_size);
    io.mapRequired("insert size deviation", data.insert_size_deviation);
    io.mapRequired("insert size left quantile", data.insert_size_left_quantile);
    io.mapRequired("insert size right quantile", data.insert_size_right_quantile);
    io.mapRequired("insert size median", data.median_insert_size);
    io.mapRequired("insert size mad", data.insert_size_mad);
    io.mapRequired("insert size distribution", data.insert_size_distribution);
    io.mapRequired("pi threshold", data.pi_threshold);
    io.mapRequired("binary reads info", data.binary_reads_info);
    io.mapRequired("single reads mapped", data.single_reads_mapped);
    io.mapRequired("library index", data.lib_index);
    io.mapRequired("number of reads", data.read_count);
    io.mapRequired("total nucleotides", data.total_nucls);
}

} }
