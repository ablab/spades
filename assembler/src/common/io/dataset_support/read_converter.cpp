//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "read_converter.hpp"

#include "assembly_graph/core/graph.hpp"
#include "io/reads/binary_streams.hpp"
#include "io/reads/rc_reader_wrapper.hpp"
#include "io/reads/multifile_reader.hpp"
#include "io/reads/converting_reader_wrapper.hpp"
#include "io/reads/edge_sequences_reader.hpp"

#include "utils/filesystem/file_opener.hpp"
#include "utils/logger/logger.hpp"

#include "threadpool/threadpool.hpp"

#include <fstream>


namespace io {

bool ReadConverter::CheckBinaryReadsExist(SequencingLibraryT& lib) {
    return fs::FileExists(lib.data().binary_reads_info.bin_reads_info_file);
}

//todo change to yaml
bool ReadConverter::LoadLibIfExists(SequencingLibraryT& lib) {
    auto& data = lib.data();

    if (!CheckBinaryReadsExist(lib))
        return false;

    auto info = fs::open_file(data.binary_reads_info.bin_reads_info_file, std::ios_base::in);
    DEBUG("Reading binary information file " << data.binary_reads_info.bin_reads_info_file);

    size_t format = 0;
    size_t lib_index = 0;

    info >> format;
    
    if (!info.eof())
        info >> lib_index;

    if (format != BINARY_FORMAT_VERSION ||
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

void ReadConverter::ConvertToBinary(SequencingLibraryT& lib,
                                    ThreadPool::ThreadPool *pool) {
    auto& data = lib.data();
    std::ofstream info;
    info.open(data.binary_reads_info.bin_reads_info_file, std::ios_base::out);
    info << "0 0 0";
    info.close();

    INFO("Converting reads to binary format for library #" << data.lib_index << " (takes a while)");
    INFO("Converting paired reads");
    BinaryWriter paired_converter(data.binary_reads_info.paired_read_prefix);

    FileReadFlags flags{ PhredOffset, /* use name */ false, /* use quality */ false, /* validate */ false };
    PairedStream paired_reader = paired_easy_reader(lib, false, 0, false, flags, pool);
    ReadStreamStat read_stat = paired_converter.ToBinary(paired_reader, lib.orientation(), pool);
    read_stat.read_count *= 2;

    INFO("Converting single reads");
    BinaryWriter single_converter(data.binary_reads_info.single_read_prefix);
    SingleStream single_reader = single_easy_reader(lib, false, false, true, flags, pool);
    read_stat.merge(single_converter.ToBinary(single_reader, pool));

    data.unmerged_read_length = read_stat.max_len;
    INFO("Converting merged reads");
    BinaryWriter merged_converter(data.binary_reads_info.merged_read_prefix);
    SingleStream merged_reader = merged_easy_reader(lib, false, true, flags, pool);
    auto merged_stats = merged_converter.ToBinary(merged_reader, pool);

    data.merged_read_length = merged_stats.max_len;
    read_stat.merge(merged_stats);
    data.read_count = read_stat.read_count;
    data.total_nucls = read_stat.total_len;

    WriteBinaryInfo(data.binary_reads_info.bin_reads_info_file, data);
}

void ReadConverter::ConvertEdgeSequencesToBinary(const debruijn_graph::Graph &g,
                                                 const std::string &contigs_output_dir, unsigned nthreads) {
    INFO("Outputting contigs to " << contigs_output_dir);

    std::unique_ptr<ThreadPool::ThreadPool> pool;
    if (nthreads > 1)
        pool = std::make_unique<ThreadPool::ThreadPool>(nthreads);

    io::BinaryWriter single_converter(fs::append_path(contigs_output_dir, "contigs"));
    io::ReadStream<io::SingleReadSeq> single_reader = io::EdgeSequencesStream(g);
    ReadStreamStat read_stat = single_converter.ToBinary(single_reader, pool.get());

    LibraryData data;
    data.lib_index = size_t(-1);
    data.unmerged_read_length = read_stat.max_len;
    data.merged_read_length = 0;
    data.read_count = read_stat.read_count;
    data.total_nucls = read_stat.total_len;
    WriteBinaryInfo(fs::append_path(contigs_output_dir, "contigs_info"), data);
}

void ReadConverter::WriteBinaryInfo(const std::string &filename, LibraryData &data) {
    std::ofstream info;
    info.open(filename, std::ios_base::out);
    info << BINARY_FORMAT_VERSION << " " <<
         data.lib_index << " " <<
         data.unmerged_read_length << " " <<
         data.merged_read_length << " " <<
         data.read_count << " " <<
         data.total_nucls << "\n";

    info.close();
    data.binary_reads_info.binary_converted = true;
}

void ConvertIfNeeded(DataSet<LibraryData> &data, unsigned nthreads) {
    std::unique_ptr<ThreadPool::ThreadPool> pool;

    if (nthreads > 1)
        pool = std::make_unique<ThreadPool::ThreadPool>(nthreads);

    for (auto &lib : data) {
        if (!ReadConverter::LoadLibIfExists(lib))
            ReadConverter::ConvertToBinary(lib, pool.get());
    }
}

BinaryPairedStreams paired_binary_readers(SequencingLibraryT &lib,
                                          bool followed_by_rc,
                                          size_t insert_size,
                                          bool include_merged) {
    const auto& data = lib.data();
    CHECK_FATAL_ERROR(data.binary_reads_info.binary_converted,
            "Lib was not converted to binary, cannot produce binary stream");

    ReadStreamList<PairedReadSeq> paired_streams;
    const size_t n = data.binary_reads_info.chunk_num;

    for (size_t i = 0; i < n; ++i) {
        ReadStream<PairedReadSeq> stream{BinaryFilePairedStream(data.binary_reads_info.paired_read_prefix,
                                                                   insert_size, n, i)};
        if (include_merged) {
            VERIFY(lib.data().unmerged_read_length != 0);
            stream = MultifileWrap<PairedReadSeq>(std::move(stream),
                                                  BinaryUnmergingPairedStream(data.binary_reads_info.merged_read_prefix,
                                                                              insert_size, lib.data().unmerged_read_length,
                                                                              n, i));
        }

        paired_streams.push_back(std::move(stream));
    }

    if (followed_by_rc)
        paired_streams = RCWrap<PairedReadSeq>(std::move(paired_streams));

    return paired_streams;
}

BinarySingleStreams single_binary_readers(SequencingLibraryT &lib,
                                          bool followed_by_rc,
                                          bool including_paired_and_merged) {
    const auto& data = lib.data();
    CHECK_FATAL_ERROR(data.binary_reads_info.binary_converted,
               "Lib was not converted to binary, cannot produce binary stream");

    BinarySingleStreams single_streams;
    const size_t n = data.binary_reads_info.chunk_num;

    for (size_t i = 0; i < n; ++i)
        single_streams.push_back(BinaryFileSingleStream(data.binary_reads_info.single_read_prefix,
                                                        n, i));

    if (including_paired_and_merged) {
        BinarySingleStreams merged_streams;
        for (size_t i = 0; i < n; ++i)
            merged_streams.push_back(BinaryFileSingleStream(data.binary_reads_info.merged_read_prefix,
                                                            n, i));
        single_streams = WrapPairsInMultifiles<SingleReadSeq>(std::move(single_streams), std::move(merged_streams));

        BinaryPairedStreams paired_streams;
        for (size_t i = 0; i < n; ++i)
            paired_streams.push_back(BinaryFilePairedStream(data.binary_reads_info.paired_read_prefix,
                                                            0, n, i));
        single_streams = WrapPairsInMultifiles<SingleReadSeq>(std::move(single_streams),
                                                              SquashingWrap<PairedReadSeq>(std::move(paired_streams)));
    }

    if (followed_by_rc)
        single_streams = RCWrap<SingleReadSeq>(std::move(single_streams));

    return single_streams;
}

BinarySingleStreams
single_binary_readers_for_libs(DataSet<LibraryData>& dataset_info,
                               const std::vector<size_t>& libs,
                               bool followed_by_rc,
                               bool including_paired_reads) {
    VERIFY(!libs.empty())
    size_t chunk_num = dataset_info[libs.front()].data().binary_reads_info.chunk_num;

    std::vector<BinarySingleStreams> streams(chunk_num);
    for (size_t i = 0; i < libs.size(); ++i) {
        VERIFY_MSG(chunk_num == dataset_info[libs[i]].data().binary_reads_info.chunk_num,
                   "Cannot create stream for multiple libraries with different chunk_num")
        BinarySingleStreams lib_streams = single_binary_readers(dataset_info[libs[i]],
                                                                followed_by_rc, including_paired_reads);

        for (size_t j = 0; j < chunk_num; ++j)
            streams[j].push_back(std::move(lib_streams[j]));
    }

    BinarySingleStreams joint_streams;
    for (size_t j = 0; j < chunk_num; ++j) {
        joint_streams.push_back(MultifileWrap<SingleReadSeq>(std::move(streams[j])));
    }
    return joint_streams;
}

}
