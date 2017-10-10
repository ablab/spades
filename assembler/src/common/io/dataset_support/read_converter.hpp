//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * read_converter.hpp
 *
 *  Created on: Apr 13, 2012
 *      Author: andrey
 */

#pragma once

#include "io/reads/binary_converter.hpp"
#include "io/reads/io_helper.hpp"
#include "dataset_readers.hpp"
#include "utils/stl_utils.hpp"

#include <fstream>

namespace io {

typedef debruijn_graph::config::dataset dataset;
typedef debruijn_graph::config::DataSetData DataSetData;
typedef SequencingLibrary<DataSetData> SequencingLibraryT;

class ReadConverter {

private:
    const static size_t current_binary_format_version = 11;

    static bool CheckBinaryReadsExist(SequencingLibraryT& lib) {
        return fs::FileExists(lib.data().binary_reads_info.bin_reads_info_file);
    }

    //FIXME change to yaml
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
            format != current_binary_format_version ||
            lib_index != data.lib_index) {
            return false;
        }

        INFO("Binary reads detected");
        info >> data.read_length;
        info >> data.merged_length;
        //info >> data.read_count;
        info >> data.total_nucls;
        data.binary_reads_info.binary_coverted = true;

        info.close();
        return true;
    }

    //FIXME change to yaml
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
        read_stat.read_count_ *= 2;

        INFO("Converting single reads");
        BinaryWriter single_converter(data.binary_reads_info.single_read_prefix,
                                          data.binary_reads_info.chunk_num,
                                          data.binary_reads_info.buffer_size);
        SingleStreamPtr single_reader = single_easy_reader(lib, false, false);
        read_stat.merge(single_converter.ToBinary(*single_reader));

        INFO("Converting merged reads");
        BinaryWriter merged_converter(data.binary_reads_info.merged_read_prefix,
                                      data.binary_reads_info.chunk_num,
                                      data.binary_reads_info.buffer_size);
        SingleStreamPtr merged_reader = merged_easy_reader(lib, false);
        auto merged_stats = merged_converter.ToBinary(*merged_reader);

        //FIXME should we put merged only stats here (or max with unmerged)
        data.merged_length = merged_stats.max_len_;
        data.read_length = read_stat.max_len_;

        //delayed merging hack to keep length stats separate
        read_stat.merge(merged_stats);
        data.read_count = read_stat.read_count_;
        data.total_nucls = read_stat.total_len_;

        //FIXME use yaml
        info.open(data.binary_reads_info.bin_reads_info_file.c_str(), std::ios_base::out);
        info << current_binary_format_version << " " <<
            data.binary_reads_info.chunk_num << " " <<
            data.lib_index << " " <<
            data.read_length << " " <<
            //data.read_count << " " <<
            data.merged_length << " " <<
            data.total_nucls << "\n";

        info.close();
        data.binary_reads_info.binary_coverted = true;
    }

public:
    static void ConvertToBinaryIfNeeded(SequencingLibraryT& lib) {
        if (lib.data().binary_reads_info.binary_coverted && CheckBinaryReadsExist(lib))
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
                                          bool include_merged,
                                          size_t read_length) {
    ReadConverter::ConvertToBinaryIfNeeded(lib);
    const auto& data = lib.data();
    VERIFY_MSG(data.binary_reads_info.binary_coverted, "Lib was not converted to binary, cannot produce binary stream");

    ReadStreamList<PairedReadSeq> paired_streams;
    const size_t n = data.binary_reads_info.chunk_num;

    for (size_t i = 0; i < n; ++i) {
        BinaryPairedStreamPtr stream = make_shared<BinaryFilePairedStream>(data.binary_reads_info.paired_read_prefix,
                                                     i, insert_size);
        if (include_merged)
            stream = MultifileWrap<PairedReadSeq>(stream,
                make_shared<BinaryUnMergedPairedStream>(data.binary_reads_info.merged_read_prefix,
                                                         i, insert_size, read_length));

        paired_streams.push_back(stream);
    }

    if (followed_by_rc) {
        paired_streams = RCWrap<PairedReadSeq>(paired_streams);
    }

    return paired_streams;
}

//FIXME discuss if we need additional flag to include merged reads
inline
BinarySingleStreams single_binary_readers(SequencingLibraryT &lib,
                                          bool followed_by_rc,
                                          bool including_paired_reads) {
    const auto& data = lib.data();
    ReadConverter::ConvertToBinaryIfNeeded(lib);
    VERIFY_MSG(data.binary_reads_info.binary_coverted, "Lib was not converted to binary, cannot produce binary stream");

    BinarySingleStreams single_streams;
    const size_t n = data.binary_reads_info.chunk_num;

    for (size_t i = 0; i < n; ++i) {
        single_streams.push_back(make_shared<BinaryFileSingleStream>(data.binary_reads_info.single_read_prefix, i));
    }

    if (including_paired_reads) {
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
        single_streams = WrapPairsInMultifiles<SingleReadSeq>(single_streams, SquashingWrap<PairedReadSeq>(paired_streams));
    }

    if (followed_by_rc) {
        single_streams = RCWrap<SingleReadSeq>(single_streams);
    }

    return single_streams;
}

//inline
////todo simplify
//BinaryPairedStreams paired_binary_readers_for_libs(dataset& dataset_info,
//                                                       const std::vector<size_t>& libs,
//                                                       bool followed_by_rc,
//                                                       size_t insert_size = 0) {
//
//    VERIFY(!libs.empty())
//    size_t chunk_num = dataset_info.reads[libs.front()].data().binary_reads_info.chunk_num;
//
//    std::vector<BinaryPairedStreams> streams(chunk_num);
//    for (size_t i = 0; i < libs.size(); ++i) {
//        VERIFY_MSG(chunk_num == dataset_info.reads[libs[i]].data().binary_reads_info.chunk_num,
//                   "Cannot create stream for multiple libraries with different chunk_num")
//        BinaryPairedStreams lib_streams = paired_binary_readers(dataset_info.reads[libs[i]], followed_by_rc, insert_size);
//        for (size_t j = 0; j < chunk_num; ++j) {
//            streams[j].push_back(lib_streams.ptr_at(j));
//        }
//    }
//
//    BinaryPairedStreams joint_streams;
//    for (size_t j = 0; j < chunk_num; ++j) {
//      joint_streams.push_back(MultifileWrap<PairedReadSeq>(streams[j]));
//    }
//    return joint_streams;
//}

inline
BinarySingleStreams single_binary_readers_for_libs(dataset& dataset_info,
                                                       const std::vector<size_t>& libs,
                                                       bool followed_by_rc,
                                                       bool including_paired_reads) {
    VERIFY(!libs.empty())
    size_t chunk_num = dataset_info.reads[libs.front()].data().binary_reads_info.chunk_num;

    std::vector<BinarySingleStreams> streams(chunk_num);
    for (size_t i = 0; i < libs.size(); ++i) {
        VERIFY_MSG(chunk_num == dataset_info.reads[libs[i]].data().binary_reads_info.chunk_num,
                   "Cannot create stream for multiple libraries with different chunk_num")
        BinarySingleStreams lib_streams = single_binary_readers(dataset_info.reads[libs[i]], followed_by_rc, including_paired_reads);

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

//inline
//BinaryPairedStreams paired_binary_readers(dataset& dataset_info,
//                                              bool followed_by_rc,
//                                              size_t insert_size = 0) {
//
//  std::vector<size_t> all_libs(dataset_info.reads.lib_count());
//  for (size_t i = 0; i < dataset_info.reads.lib_count(); ++i) {
//      all_libs[i] = i;
//  }
//  return paired_binary_readers_for_libs(dataset_info, all_libs, followed_by_rc, insert_size);
//}
//
//inline
//BinarySingleStreams single_binary_readers(dataset& dataset_info,
//                                              bool followed_by_rc,
//                                              bool including_paired_reads) {
//  std::vector<size_t> all_libs(dataset_info.reads.lib_count());
//  for (size_t i = 0; i < dataset_info.reads.lib_count(); ++i) {
//    all_libs[i] = i;
//  }
//  return single_binary_readers_for_libs(dataset_info, all_libs, followed_by_rc, including_paired_reads);
//}
//
//inline
//BinarySingleStreamPtr single_binary_multireader(dataset& dataset_info, bool followed_by_rc, bool including_paired_reads) {
//    return MultifileWrap<SingleReadSeq>(single_binary_readers(dataset_info, followed_by_rc, including_paired_reads));
//}
//
//inline
//BinaryPairedStreamPtr paired_binary_multireader(dataset& dataset_info, bool followed_by_rc, size_t insert_size = 0) {
//    return MultifileWrap<PairedReadSeq>(paired_binary_readers(dataset_info, followed_by_rc, insert_size));
//}


}
