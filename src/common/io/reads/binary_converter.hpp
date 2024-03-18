//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "read_stream.hpp"
#include "single_read.hpp"
#include "paired_read.hpp"
#include "orientation.hpp"

#include "library/library_fwd.hpp"

#include <fstream>

namespace ThreadPool {
class ThreadPool;
};

namespace io {

template<class Read>
using ReadTagger = std::function<uint64_t(const Read&)>;

class BinaryWriter {
    const std::string file_name_prefix_;
    std::unique_ptr<std::ofstream> file_ds_, offset_ds_;

    template<class Writer, class Read>
    ReadStreamStat ToBinary(const Writer &writer, io::ReadStream<Read> &stream,
                            ThreadPool::ThreadPool *pool = nullptr);

    template<class Read>
    struct TrivialTagger {
        uint64_t operator()(const Read &) const { return 0; }
    };

    template<class Read>
    struct IdempotentTagger {
        uint64_t operator()(const Read &r) const { return r.tag(); }
    };

    template<class Read>
    struct BarcodeTagger {
        uint64_t operator()(const Read &r) const { return r.aux().sequence().data()[0]; }
    };

public:
    typedef size_t CountType;
    static constexpr size_t CHUNK = 100;
    static constexpr size_t BUF_SIZE = 50000;

    BinaryWriter(const std::string &file_name_prefix);

    ~BinaryWriter() = default;

    ReadStreamStat ToBinary(io::ReadStream<io::SingleReadSeq>& stream,
                            ThreadPool::ThreadPool *pool = nullptr,
                            ReadTagger<io::SingleReadSeq> tagger = IdempotentTagger<io::SingleReadSeq>());
    ReadStreamStat ToBinary(io::ReadStream<io::SingleRead>& stream,
                            ThreadPool::ThreadPool *pool = nullptr,
                            ReadTagger<io::SingleRead> tagger = TrivialTagger<io::SingleRead>());
    ReadStreamStat ToBinary(io::ReadStream<io::PairedRead>& stream,
                            LibraryOrientation orientation = LibraryOrientation::Undefined,
                            ThreadPool::ThreadPool *pool = nullptr,
                            ReadTagger<io::SingleRead> tagger = TrivialTagger<io::SingleRead>());
    ReadStreamStat ToBinary(io::ReadStream<io::TellSeqRead>& stream,
                            LibraryOrientation orientation = LibraryOrientation::Undefined,
                            ThreadPool::ThreadPool *pool = nullptr,
                            ReadTagger<io::TellSeqRead> tagger = BarcodeTagger<io::TellSeqRead>());
};

}
