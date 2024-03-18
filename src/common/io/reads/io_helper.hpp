//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "read_stream_vector.hpp"
#include "single_read.hpp"
#include "paired_read.hpp"

#include "library/library_fwd.hpp"


namespace ThreadPool {
class ThreadPool;
}

namespace io {
typedef ReadStream<SingleRead> SingleStream;
typedef ReadStreamList<SingleRead> SingleStreams;

typedef ReadStream<PairedRead> PairedStream;
typedef ReadStreamList<PairedRead> PairedStreams;

typedef ReadStream<TellSeqRead> TellSeqStream;
typedef ReadStreamList<TellSeqRead> TellSeqStreams;

typedef ReadStream<SingleReadSeq> BinarySingleStream;
typedef ReadStreamList<SingleReadSeq> BinarySingleStreams;

typedef ReadStream<PairedReadSeq> BinaryPairedStream;
typedef ReadStreamList<PairedReadSeq> BinaryPairedStreams;

SingleStream EasyStream(const std::filesystem::path& filename, bool followed_by_rc,
                        bool handle_Ns = true,
                        FileReadFlags flags = FileReadFlags(),
                        ThreadPool::ThreadPool *pool = nullptr);

PairedStream EasyWrapPairedStream(PairedStream stream,
                                  bool followed_by_rc,
                                  LibraryOrientation orientation, bool handle_Ns=true);

PairedStream PairedEasyStream(const std::filesystem::path& filename1, const std::filesystem::path& filename2,
                              bool followed_by_rc, size_t insert_size,
                              bool use_orientation = true, bool handle_Ns=true,
                              LibraryOrientation orientation = LibraryOrientation::FR,
                              FileReadFlags flags = FileReadFlags(),
                              ThreadPool::ThreadPool *pool = nullptr);

PairedStream PairedEasyStream(const std::filesystem::path& filename, bool followed_by_rc,
                              size_t insert_size,
                              bool use_orientation = true,  bool handle_Ns=true,
                              LibraryOrientation orientation = LibraryOrientation::FR,
                              FileReadFlags flags = FileReadFlags(),
                              ThreadPool::ThreadPool *pool = nullptr);

TellSeqStream TellSeqEasyStream(const std::filesystem::path& filename1, const std::filesystem::path& filename2,
                                const std::filesystem::path& aux,
                                bool followed_by_rc, size_t insert_size,
                                bool use_orientation = true, bool handle_Ns=true,
                                LibraryOrientation orientation = LibraryOrientation::FR,
                                FileReadFlags flags = FileReadFlags(),
                                ThreadPool::ThreadPool *pool = nullptr);

}
