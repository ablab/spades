//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/config_struct.hpp"
#include "io/reads/binary_converter.hpp"
#include "io/reads/io_helper.hpp"
#include "dataset_readers.hpp"

namespace io {

typedef debruijn_graph::config::LibraryData LibraryData;
typedef SequencingLibrary<LibraryData> SequencingLibraryT;

class ReadConverter {
    static constexpr size_t BINARY_FORMAT_VERSION = 13;

    static bool CheckBinaryReadsExist(SequencingLibraryT& lib);

public:
    static bool LoadLibIfExists(SequencingLibraryT& lib);
    static void ConvertToBinary(SequencingLibraryT& lib,
                                ThreadPool::ThreadPool *pool = nullptr);
};

void ConvertIfNeeded(DataSet<LibraryData> &data, unsigned nthreads);

BinaryPairedStreams paired_binary_readers(SequencingLibraryT &lib,
                                          bool followed_by_rc,
                                          size_t insert_size,
                                          bool include_merged);
BinarySingleStreams single_binary_readers(SequencingLibraryT &lib,
                                          bool followed_by_rc,
                                          bool including_paired_and_merged);

BinarySingleStreams single_binary_readers_for_libs(DataSet<LibraryData>& dataset_info,
                                                   const std::vector<size_t>& libs,
                                                   bool followed_by_rc = true,
                                                   bool including_paired_reads = true);

}
