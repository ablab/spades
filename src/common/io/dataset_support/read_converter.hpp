//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "dataset_readers.hpp"

#include "io/reads/single_read.hpp"
#include "io/reads/binary_converter.hpp"
#include "io/reads/io_helper.hpp"

#include "configs/config_struct.hpp"

#include <functional>

namespace debruijn_graph {
class DeBruijnGraph;
typedef DeBruijnGraph Graph;
}

namespace io {

typedef debruijn_graph::config::LibraryData LibraryData;
typedef SequencingLibrary<LibraryData> SequencingLibraryT;

class ReadConverter {
    static constexpr size_t BINARY_FORMAT_VERSION = 14;

    static bool CheckBinaryReadsExist(SequencingLibraryT& lib);
    static void WriteBinaryInfo(const std::filesystem::path& filename, LibraryData& data);
public:
    struct TrivialTagger {
        uint64_t operator()(const io::SingleRead &) const { return 0; }
    };
    
    static bool LoadLibIfExists(SequencingLibraryT& lib);
    static void ConvertToBinary(SequencingLibraryT& lib,
                                ThreadPool::ThreadPool *pool = nullptr,
                                FileReadFlags flags = FileReadFlags::empty(),
                                ReadTagger<io::SingleRead> tagger = TrivialTagger());

    static void ConvertEdgeSequencesToBinary(const debruijn_graph::Graph &g, const std::filesystem::path &contigs_output_dir,
                                             unsigned nthreads);
};

void ConvertIfNeeded(DataSet<LibraryData> &data, unsigned nthreads = 1,
                     FileReadFlags flags = FileReadFlags::empty(),
                     ReadTagger<io::SingleRead> tagger = ReadConverter::TrivialTagger());

BinaryPairedStreams paired_binary_readers(SequencingLibraryT &lib,
                                          bool followed_by_rc,
                                          size_t insert_size,
                                          bool include_merged,
                                          size_t chunk_num = 0);
BinarySingleStreams single_binary_readers(SequencingLibraryT &lib,
                                          bool followed_by_rc,
                                          bool including_paired_and_merged,
                                          size_t chunk_num = 0);

BinarySingleStreams single_binary_readers_for_libs(DataSet<LibraryData>& dataset_info,
                                                   const std::vector<size_t>& libs,
                                                   bool followed_by_rc = true,
                                                   bool including_paired_reads = true,
                                                   size_t chunk_num = 0);

}
