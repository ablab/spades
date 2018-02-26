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

namespace io {
    typedef ReadStream<SingleRead> SingleStream;
    typedef std::shared_ptr<SingleStream> SingleStreamPtr;
    typedef ReadStreamList<SingleRead> SingleStreams;

    typedef ReadStream<PairedRead> PairedStream;
    typedef std::shared_ptr<PairedStream> PairedStreamPtr;
    typedef ReadStreamList<PairedRead> PairedStreams;

    typedef ReadStream<SingleReadSeq> BinarySingleStream;
    typedef std::shared_ptr<BinarySingleStream> BinarySingleStreamPtr;
    typedef ReadStreamList<SingleReadSeq> BinarySingleStreams;

    typedef ReadStream<PairedReadSeq> BinaryPairedStream;
    typedef std::shared_ptr<BinaryPairedStream> BinaryPairedStreamPtr;
    typedef ReadStreamList<PairedReadSeq> BinaryPairedStreams;

    inline SingleStreamPtr EasyStream(const std::string& filename, bool followed_by_rc,
                                      bool handle_Ns = true, OffsetType offset_type = PhredOffset) {
        SingleStreamPtr reader = make_shared<FileReadStream>(filename, offset_type);
        if (handle_Ns) {
            reader = LongestValidWrap<SingleRead>(reader);
        }
        if (followed_by_rc) {
            reader = RCWrap<SingleRead>(reader);
        }
        return reader;
    }

    inline PairedStreamPtr EasyWrapPairedStream(PairedStreamPtr reader,
                                                bool followed_by_rc,
                                                LibraryOrientation orientation) {
        reader = make_shared<OrientationChangingWrapper<PairedRead>>(reader, orientation);
        reader = LongestValidWrap<PairedRead>(reader);
        if (followed_by_rc) {
            reader = RCWrap<PairedRead>(reader);
        }
        return reader;
    }

    inline PairedStreamPtr PairedEasyStream(const std::string& filename1, const std::string& filename2,
                                     bool followed_by_rc, size_t insert_size,
                                     bool use_orientation = true, LibraryOrientation orientation = LibraryOrientation::FR,
                                     OffsetType offset_type = PhredOffset) {
        PairedStreamPtr reader = make_shared<SeparatePairedReadStream>(filename1, filename2,
                                                                       insert_size, offset_type);

        return EasyWrapPairedStream(reader, followed_by_rc,
                                    use_orientation ? orientation : LibraryOrientation::Undefined);
    }

    inline PairedStreamPtr PairedEasyStream(const std::string& filename, bool followed_by_rc,
            size_t insert_size,
            bool use_orientation = true, LibraryOrientation orientation = LibraryOrientation::FR,
            OffsetType offset_type = PhredOffset) {
        PairedStreamPtr reader = make_shared<InterleavingPairedReadStream>(filename,
                                                                           insert_size,
                                                                           offset_type);
        return EasyWrapPairedStream(reader, followed_by_rc,
                                    use_orientation ? orientation : LibraryOrientation::Undefined);
    }
}
