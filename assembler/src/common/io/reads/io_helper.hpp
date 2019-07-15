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
#include "file_reader.hpp"
#include "paired_readers.hpp"
#include "binary_streams.hpp"
#include "multifile_reader.hpp"
#include "converting_reader_wrapper.hpp"
#include "longest_valid_wrapper.hpp"
#include "rc_reader_wrapper.hpp"
#include "async_read_stream.hpp"

namespace io {
    typedef ReadStream<SingleRead> SingleStream;
    typedef ReadStreamList<SingleRead> SingleStreams;

    typedef ReadStream<PairedRead> PairedStream;
    typedef ReadStreamList<PairedRead> PairedStreams;

    typedef ReadStream<SingleReadSeq> BinarySingleStream;
    typedef ReadStreamList<SingleReadSeq> BinarySingleStreams;

    typedef ReadStream<PairedReadSeq> BinaryPairedStream;
    typedef ReadStreamList<PairedReadSeq> BinaryPairedStreams;

    inline SingleStream EasyStream(const std::string& filename, bool followed_by_rc,
                                   bool handle_Ns = true, OffsetType offset_type = PhredOffset,
                                   ThreadPool::ThreadPool *pool = nullptr) {
        SingleStream reader  = (pool ?
                                make_async_stream<FileReadStream>(*pool, filename, offset_type) :
                                FileReadStream(filename, offset_type));
        if (handle_Ns)
            reader = LongestValidWrap<SingleRead>(std::move(reader));
        if (followed_by_rc)
            reader = RCWrap<SingleRead>(std::move(reader));

        return reader;
    }

    inline PairedStream EasyWrapPairedStream(PairedStream stream,
                                             bool followed_by_rc,
                                             LibraryOrientation orientation) {
        PairedStream reader{OrientationChangingWrapper<PairedRead>(std::move(stream), orientation)};
        reader = LongestValidWrap<PairedRead>(std::move(reader));
        if (followed_by_rc)
            reader = RCWrap<PairedRead>(std::move(reader));

        return reader;
    }

    inline PairedStream PairedEasyStream(const std::string& filename1, const std::string& filename2,
                                         bool followed_by_rc, size_t insert_size,
                                         bool use_orientation = true, LibraryOrientation orientation = LibraryOrientation::FR,
                                         OffsetType offset_type = PhredOffset,
                                         ThreadPool::ThreadPool *pool = nullptr) {
        return EasyWrapPairedStream(SeparatePairedReadStream(filename1, filename2, insert_size, offset_type,
                                                             pool),
                                    followed_by_rc,
                                    use_orientation ? orientation : LibraryOrientation::Undefined);
    }

    inline PairedStream PairedEasyStream(const std::string& filename, bool followed_by_rc,
                                         size_t insert_size,
                                         bool use_orientation = true, LibraryOrientation orientation = LibraryOrientation::FR,
                                         OffsetType offset_type = PhredOffset,
                                         ThreadPool::ThreadPool *pool = nullptr) {
        return EasyWrapPairedStream(InterleavingPairedReadStream(filename, insert_size, offset_type, pool), followed_by_rc,
                                    use_orientation ? orientation : LibraryOrientation::Undefined);
    }
}
