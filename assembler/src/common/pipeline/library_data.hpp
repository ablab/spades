//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <map>
#include <string>

// Forward decls for LLVM YAML API
namespace llvm { namespace yaml { class IO; template<typename T> struct MappingTraits; } }

namespace debruijn_graph {

namespace config {

struct LibraryData {
    size_t unmerged_read_length;
    size_t merged_read_length;
    double mean_insert_size;
    double insert_size_deviation;
    double insert_size_left_quantile;
    double insert_size_right_quantile;
    double median_insert_size;
    double insert_size_mad;
    std::map<int, size_t> insert_size_distribution;

    size_t lib_index;
    bool single_reads_mapped;
    uint64_t total_nucls;
    size_t read_count;

    double pi_threshold;

    struct BinaryReadsInfo {
        BinaryReadsInfo() {}

        bool binary_converted = false;
        std::string bin_reads_info_file;
        std::string paired_read_prefix;
        std::string merged_read_prefix;
        std::string single_read_prefix;
        size_t chunk_num = 0;
    } binary_reads_info;

    LibraryData()
            : unmerged_read_length(0), merged_read_length(0),
              mean_insert_size(0.0),
              insert_size_deviation(0.0), insert_size_left_quantile(0.0), insert_size_right_quantile(0.0),
              median_insert_size(0.0), insert_size_mad(0.0),
              lib_index(0),
              single_reads_mapped(false),
              total_nucls(0), read_count(0),
              pi_threshold(0.0) {}
};

} // namespace config

} // namespace debruijn_graph

namespace llvm { namespace yaml {

template<>
struct MappingTraits<debruijn_graph::config::LibraryData::BinaryReadsInfo> {
    static void mapping(IO &io, debruijn_graph::config::LibraryData::BinaryReadsInfo &info);
};

template<>
struct MappingTraits<debruijn_graph::config::LibraryData> {
    static void mapping(IO &io, debruijn_graph::config::LibraryData &data);
};

} }
