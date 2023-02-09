//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "io_helper.hpp"

#include "file_reader.hpp"
#include "io/reads/paired_read.hpp"
#include "paired_readers.hpp"
#include "multifile_reader.hpp"
#include "converting_reader_wrapper.hpp"
#include "longest_valid_wrapper.hpp"
#include "rc_reader_wrapper.hpp"
#include "async_read_stream.hpp"

namespace io {

SingleStream EasyStream(const std::filesystem::path& filename, bool followed_by_rc,
                        bool handle_Ns,
                        FileReadFlags flags,
                        ThreadPool::ThreadPool *pool) {
    SingleStream reader  = (pool ?
                            make_async_stream<FileReadStream>(*pool, filename, flags) :
                            FileReadStream(filename, flags));
    if (handle_Ns)
        reader = LongestValidWrap<SingleRead>(std::move(reader));
    if (followed_by_rc)
        reader = RCWrap<SingleRead>(std::move(reader));

    return reader;
}

PairedStream EasyWrapPairedStream(PairedStream stream,
                                  bool followed_by_rc,
                                  LibraryOrientation orientation,
                                  bool handle_Ns) {
    PairedStream reader{std::move(stream)};
    if (orientation != LibraryOrientation::Undefined)
        reader = OrientationChangingWrapper<PairedRead>(std::move(reader), orientation);
    if (handle_Ns) {
        reader = LongestValidWrap<PairedRead>(std::move(reader));
    }
    if (followed_by_rc)
        reader = RCWrap<PairedRead>(std::move(reader));

    return reader;
}

TellSeqStream EasyWrapTellSeqStream(TellSeqStream stream,
                                    bool followed_by_rc,
                                    LibraryOrientation orientation,
                                    bool handle_Ns) {
    TellSeqStream reader{std::move(stream)};
    if (orientation != LibraryOrientation::Undefined)
        reader = OrientationChangingWrapper<TellSeqRead>(std::move(reader), orientation);
    if (handle_Ns) {
        reader = LongestValidWrap<TellSeqRead>(std::move(reader));
    }
    if (followed_by_rc)
        reader = RCWrap<TellSeqRead>(std::move(reader));

    return reader;
}

PairedStream PairedEasyStream(const std::filesystem::path& filename1, const std::filesystem::path& filename2,
                              bool followed_by_rc, size_t insert_size,
                              bool use_orientation, bool handle_Ns, LibraryOrientation orientation,
                              FileReadFlags flags,
                              ThreadPool::ThreadPool *pool) {
    return EasyWrapPairedStream(SeparatePairedReadStream(filename1, filename2, insert_size, flags,
                                                         pool),
                                followed_by_rc,
                                use_orientation ? orientation : LibraryOrientation::Undefined, handle_Ns);
}


TellSeqStream TellSeqEasyStream(const std::filesystem::path& filename1, const std::filesystem::path& filename2,
                                const std::filesystem::path& aux,
                                bool followed_by_rc, size_t insert_size,
                                bool use_orientation, bool handle_Ns, LibraryOrientation orientation,
                                FileReadFlags flags,
                                ThreadPool::ThreadPool *pool) {
    TellSeqReadStream s(filename1, filename2, aux,
                        insert_size, flags,  pool);

    return EasyWrapTellSeqStream(TellSeqReadStream(filename1, filename2, aux,
                                                   insert_size, flags, pool),
                                followed_by_rc,
                                use_orientation ? orientation : LibraryOrientation::Undefined, handle_Ns);
}

PairedStream PairedEasyStream(const std::filesystem::path& filename, bool followed_by_rc,
                              size_t insert_size,
                              bool use_orientation, bool handle_Ns, LibraryOrientation orientation,
                              FileReadFlags flags,
                              ThreadPool::ThreadPool *pool) {
    return EasyWrapPairedStream(InterleavingPairedReadStream(filename, insert_size, flags, pool), followed_by_rc,
                                use_orientation ? orientation : LibraryOrientation::Undefined, handle_Ns);
}

}
