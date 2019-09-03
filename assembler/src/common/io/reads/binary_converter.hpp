//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "read_stream.hpp"
#include "single_read.hpp"
#include "paired_read.hpp"
#include "orientation.hpp"

#include "pipeline/library_fwd.hpp"

#include <fstream>

namespace ThreadPool {
class ThreadPool;
};

namespace io {

class BinaryWriter {
    const std::string file_name_prefix_;
    std::unique_ptr<std::ofstream> file_ds_, offset_ds_;

    template<class Writer, class Read>
    ReadStreamStat ToBinary(const Writer &writer, io::ReadStream<Read> &stream,
                            ThreadPool::ThreadPool *pool = nullptr);

public:
    typedef size_t CountType;
    static constexpr size_t CHUNK = 100;
    static constexpr size_t BUF_SIZE = 50000;

    BinaryWriter(const std::string &file_name_prefix);

    ~BinaryWriter() = default;

    ReadStreamStat ToBinary(io::ReadStream<io::SingleReadSeq>& stream,
                            ThreadPool::ThreadPool *pool = nullptr);
    ReadStreamStat ToBinary(io::ReadStream<io::SingleRead>& stream,
                            ThreadPool::ThreadPool *pool = nullptr);
    ReadStreamStat ToBinary(io::ReadStream<io::PairedReadSeq>& stream,
                            LibraryOrientation orientation = LibraryOrientation::Undefined,
                            ThreadPool::ThreadPool *pool = nullptr);
    ReadStreamStat ToBinary(io::ReadStream<io::PairedRead>& stream,
                            LibraryOrientation orientation = LibraryOrientation::Undefined,
                            ThreadPool::ThreadPool *pool = nullptr);
};

}
