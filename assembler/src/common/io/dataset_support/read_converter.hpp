//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io/reads/binary_converter.hpp"
#include "io/reads/io_helper.hpp"
#include "dataset_readers.hpp"
#include "utils/stl_utils.hpp"

#include <fstream>

namespace io {

typedef debruijn_graph::config::LibraryData LibraryData;
typedef SequencingLibrary<LibraryData> SequencingLibraryT;

class ReadConverter {

    static const size_t BINARY_FORMAT_VERSION = 12;

    static bool CheckBinaryReadsExist(SequencingLibraryT& lib) {
        return fs::FileExists(lib.data().binary_reads_info.bin_reads_info_file);
    }

    //todo change to yaml
    static bool LoadLibIfExists(SequencingLibraryT& lib) {
        auto& data = lib.data();

        if (!CheckBinaryReadsExist(lib))
            return false;

        std::ifstream info;
        info.open(data.binary_reads_info.bin_reads_info_file.c_str(), std::ios_base::in);
        DEBUG("Reading binary information file " << data.binary_reads_info.bin_reads_info_file);

        size_t chunk_num = 0;
        size_t format = 0;
        size_t lib_index = 0;

        info >> format;
        if (!info.eof()) {
            info >> chunk_num;
        }
        if (!info.eof()) {
            info >> lib_index;
        }

        if (chunk_num != data.binary_reads_info.chunk_num ||
            format != BINARY_FORMAT_VERSION ||
            lib_index != data.lib_index) {
            return false;
        }

        INFO("Binary reads detected");
        info >> data.unmerged_read_length;
        info >> data.merged_read_length;
        info >> data.read_count;
        info >> data.total_nucls;
        data.binary_reads_info.binary_converted = true;

        info.close();
        return true;
    }

    static void ConvertToBinary(SequencingLibraryT& lib) {
        auto& data = lib.data();
        std::ofstream info;
        info.open(data.binary_reads_info.bin_reads_info_file.c_str(), std::ios_base::out);
        info << "0 0 0";
        info.close();

        INFO("Converting reads to binary format for library #" << data.lib_index << " (takes a while)");
        INFO("Converting paired reads");
        BinaryWriter paired_converter(data.binary_reads_info.paired_read_prefix,
                                          data.binary_reads_info.chunk_num,
                                          data.binary_reads_info.buffer_size);

        PairedStreamPtr paired_reader = paired_easy_reader(lib, false, 0, false, PhredOffset);
        ReadStreamStat read_stat = paired_converter.ToBinary(*paired_reader, lib.orientation());
        read_stat.read_count *= 2;

        INFO("Converting single reads");
        BinaryWriter single_converter(data.binary_reads_info.single_read_prefix,
                                          data.binary_reads_info.chunk_num,
                                          data.binary_reads_info.buffer_size);
        SingleStreamPtr single_reader = single_easy_reader(lib, false, false);
        read_stat.merge(single_converter.ToBinary(*single_reader));

        data.unmerged_read_length = read_stat.max_len;
        INFO("Converting merged reads");
        BinaryWriter merged_converter(data.binary_reads_info.merged_read_prefix,
                                      data.binary_reads_info.chunk_num,
                                      data.binary_reads_info.buffer_size);
        SingleStreamPtr merged_reader = merged_easy_reader(lib, false);
        auto merged_stats = merged_converter.ToBinary(*merged_reader);

        data.merged_read_length = merged_stats.max_len;
        read_stat.merge(merged_stats);
        data.read_count = read_stat.read_count;
        data.total_nucls = read_stat.total_len;

        //todo use yaml or reuse same data!!!
        info.open(data.binary_reads_info.bin_reads_info_file.c_str(), std::ios_base::out);
        info << BINARY_FORMAT_VERSION << " " <<
            data.binary_reads_info.chunk_num << " " <<
            data.lib_index << " " <<
            data.unmerged_read_length << " " <<
            data.merged_read_length << " " <<
            data.read_count << " " <<
            data.total_nucls << "\n";

        info.close();
        data.binary_reads_info.binary_converted = true;
    }

public:
    static void ConvertToBinaryIfNeeded(SequencingLibraryT& lib) {
        if (lib.data().binary_reads_info.binary_converted && CheckBinaryReadsExist(lib))
            return;

        if (LoadLibIfExists(lib)) {
            return;
        }

        ConvertToBinary(lib);
    }
};

inline
BinaryPairedStreams paired_binary_readers(SequencingLibraryT &lib,
                                          bool followed_by_rc,
                                          size_t insert_size,
                                          bool include_merged) {
    ReadConverter::ConvertToBinaryIfNeeded(lib);
    const auto& data = lib.data();
    CHECK_FATAL_ERROR(data.binary_reads_info.binary_converted, 
            "Lib was not converted to binary, cannot produce binary stream");

    ReadStreamList<PairedReadSeq> paired_streams;
    const size_t n = data.binary_reads_info.chunk_num;

    for (size_t i = 0; i < n; ++i) {
        BinaryPairedStreamPtr stream = make_shared<BinaryFilePairedStream>(data.binary_reads_info.paired_read_prefix,
                                                     i, insert_size);
        if (include_merged) {
            VERIFY(lib.data().unmerged_read_length != 0);
            stream = MultifileWrap<PairedReadSeq>(stream,
                                                  make_shared<BinaryUnmergingPairedStream>(data.binary_reads_info.merged_read_prefix,
                                                                                           i, insert_size, lib.data().unmerged_read_length));
        }

        paired_streams.push_back(stream);
    }

    if (followed_by_rc) {
        paired_streams = RCWrap<PairedReadSeq>(paired_streams);
    }

    return paired_streams;
}

inline
BinarySingleStreams single_binary_readers(SequencingLibraryT &lib,
                                          bool followed_by_rc,
                                          bool including_paired_and_merged) {
    const auto& data = lib.data();
    ReadConverter::ConvertToBinaryIfNeeded(lib);
    CHECK_FATAL_ERROR(data.binary_reads_info.binary_converted,
               "Lib was not converted to binary, cannot produce binary stream");

    BinarySingleStreams single_streams;
    const size_t n = data.binary_reads_info.chunk_num;

    for (size_t i = 0; i < n; ++i) {
        single_streams.push_back(make_shared<BinaryFileSingleStream>(data.binary_reads_info.single_read_prefix, i));
    }

    if (including_paired_and_merged) {
        BinarySingleStreams merged_streams;
        for (size_t i = 0; i < n; ++i) {
            merged_streams.push_back(make_shared<BinaryFileSingleStream>(data.binary_reads_info.merged_read_prefix, i));
        }
        single_streams = WrapPairsInMultifiles<SingleReadSeq>(single_streams, merged_streams);

        BinaryPairedStreams paired_streams;
        for (size_t i = 0; i < n; ++i) {
            paired_streams.push_back(make_shared<BinaryFilePairedStream>(data.binary_reads_info.paired_read_prefix,
                                                                             i, 0));
        }
        single_streams = WrapPairsInMultifiles<SingleReadSeq>(single_streams,
                                                              SquashingWrap<PairedReadSeq>(paired_streams));
    }

    if (followed_by_rc) {
        single_streams = RCWrap<SingleReadSeq>(single_streams);
    }

    return single_streams;
}

inline BinarySingleStreams
single_binary_readers_for_libs(DataSet<LibraryData>& dataset_info,
                               const std::vector<size_t>& libs,
                               bool followed_by_rc = true,
                               bool including_paired_reads = true) {
    VERIFY(!libs.empty())
    size_t chunk_num = dataset_info[libs.front()].data().binary_reads_info.chunk_num;

    std::vector<BinarySingleStreams> streams(chunk_num);
    for (size_t i = 0; i < libs.size(); ++i) {
        VERIFY_MSG(chunk_num == dataset_info[libs[i]].data().binary_reads_info.chunk_num,
                   "Cannot create stream for multiple libraries with different chunk_num")
        BinarySingleStreams lib_streams = single_binary_readers(dataset_info[libs[i]],
                                                                followed_by_rc, including_paired_reads);

        for (size_t j = 0; j < chunk_num; ++j) {
            streams[j].push_back(lib_streams.ptr_at(j));
        }
    }

    BinarySingleStreams joint_streams;
    for (size_t j = 0; j < chunk_num; ++j) {
        joint_streams.push_back(MultifileWrap<SingleReadSeq>(streams[j]));
    }
    return joint_streams;
}

}
