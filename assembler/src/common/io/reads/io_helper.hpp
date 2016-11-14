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
#include "careful_filtering_reader_wrapper.hpp"
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

    inline BinarySingleStreams apply_single_wrappers(bool followed_by_rc,
                                                     BinarySingleStreams& single_readers,
                                                     BinaryPairedStreams* paired_readers = 0) {
        VERIFY(single_readers.size() != 0);
        BinarySingleStreams readers = single_readers;

        if (paired_readers != 0) {
            VERIFY(single_readers.size() == paired_readers->size());
            BinarySingleStreams squashed_paired = SquashingWrap<PairedReadSeq>(*paired_readers);
            readers = WrapPairsInMultifiles<SingleReadSeq>(squashed_paired, readers);
        }

        if (followed_by_rc) {
            readers = RCWrap<SingleReadSeq>(readers);
        }
        return readers;
    }

    //todo make deprecated
    inline BinaryPairedStreams apply_paired_wrappers(bool followed_by_rc,
                                                     BinaryPairedStreams& readers) {
        VERIFY(readers.size() != 0);
        if (followed_by_rc) {
            return RCWrap<PairedReadSeq>(readers);
        } else {
            return readers;
        }
    }
    
    inline SingleStreamPtr EasyStream(const std::string& filename, bool followed_by_rc,
                                      bool handle_Ns = true, OffsetType offset_type = PhredOffset) {
        SingleStreamPtr reader = make_shared<FileReadStream>(filename, offset_type);
        if (handle_Ns) {
            reader = CarefulFilteringWrap<SingleRead>(reader);
        }
        if (followed_by_rc) {
            reader = RCWrap<SingleRead>(reader);
        }
        return reader;
    }

    inline PairedStreamPtr WrapPairedStream(PairedStreamPtr reader,
                                            bool followed_by_rc,
                                            bool use_orientation = false,
                                            LibraryOrientation orientation = LibraryOrientation::Undefined) {
        PairedStreamPtr answer = reader;
        answer = CarefulFilteringWrap<PairedRead>(answer, use_orientation, orientation);
        if (followed_by_rc) {
            answer = RCWrap<PairedRead>(answer);
        }
        return answer;

    }

    inline PairedStreamPtr PairedEasyStream(const std::string& filename1, const std::string& filename2,
                                     bool followed_by_rc, size_t insert_size, bool change_read_order = false,
                                     bool use_orientation = true, LibraryOrientation orientation = LibraryOrientation::FR,
                                     OffsetType offset_type = PhredOffset) {
        PairedStreamPtr reader = make_shared<SeparatePairedReadStream>(filename1, filename2, insert_size,
                                                             change_read_order, use_orientation,
                                                             orientation, offset_type);
        //Use orientation for IS calculation if it's not done by changer
        return WrapPairedStream(reader, followed_by_rc, !use_orientation, orientation);
    }

    inline PairedStreamPtr PairedEasyStream(const std::string& filename, bool followed_by_rc,
            size_t insert_size, bool change_read_order = false,
            bool use_orientation = true, LibraryOrientation orientation = LibraryOrientation::FR,
            OffsetType offset_type = PhredOffset) {
        PairedStreamPtr reader = make_shared<InterleavingPairedReadStream>(filename, insert_size, change_read_order,
                                use_orientation, orientation, offset_type);
        //Use orientation for IS calculation if it's not done by changer
        return WrapPairedStream(reader, followed_by_rc, !use_orientation, orientation);
    }
}
