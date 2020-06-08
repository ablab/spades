//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "read_stream_vector.hpp"
#include "single_read.hpp"
#include "paired_read.hpp"

#include "pipeline/library_fwd.hpp"

#include <string>

namespace ThreadPool {
class ThreadPool;
}

namespace io {
typedef ReadStream<SingleRead> SingleStream;
typedef ReadStreamList<SingleRead> SingleStreams;

typedef ReadStream<PairedRead> PairedStream;
typedef ReadStreamList<PairedRead> PairedStreams;

typedef ReadStream<SingleReadSeq> BinarySingleStream;
typedef ReadStreamList<SingleReadSeq> BinarySingleStreams;

typedef ReadStream<PairedReadSeq> BinaryPairedStream;
typedef ReadStreamList<PairedReadSeq> BinaryPairedStreams;

SingleStream EasyStream(const std::string& filename, bool followed_by_rc,
                        bool handle_Ns = true,
                        FileReadFlags flags = FileReadFlags(),
                        ThreadPool::ThreadPool *pool = nullptr);
PairedStream EasyWrapPairedStream(PairedStream stream,
                                  bool followed_by_rc,
                                  LibraryOrientation orientation);
PairedStream PairedEasyStream(const std::string& filename1, const std::string& filename2,
                              bool followed_by_rc, size_t insert_size,
                              bool use_orientation = true, LibraryOrientation orientation = LibraryOrientation::FR,
                              FileReadFlags flags = FileReadFlags(),
                              ThreadPool::ThreadPool *pool = nullptr);
PairedStream PairedEasyStream(const std::string& filename, bool followed_by_rc,
                              size_t insert_size,
                              bool use_orientation = true, LibraryOrientation orientation = LibraryOrientation::FR,
                              FileReadFlags flags = FileReadFlags(),
                              ThreadPool::ThreadPool *pool = nullptr);
}
