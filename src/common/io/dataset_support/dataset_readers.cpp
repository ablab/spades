//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2019-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************


#include "dataset_readers.hpp"

#include "io/reads/file_read_flags.hpp"
#include "io/reads/io_helper.hpp"
#include "io/reads/multifile_reader.hpp"

#include "io/reads/paired_read.hpp"
#include "library/library.hpp"

namespace io {

PairedStream paired_easy_reader(const SequencingLibraryBase &lib,
                                bool followed_by_rc,
                                size_t insert_size,
                                bool use_orientation,
                                bool handle_Ns,
                                FileReadFlags flags,
                                ThreadPool::ThreadPool *pool) {
    return MultifileWrap<PairedRead>(paired_easy_readers(lib, followed_by_rc, insert_size,
                                                         use_orientation, handle_Ns, flags, pool));
}


TellSeqStream tellseq_easy_reader(const SequencingLibraryBase &lib,
                                  bool followed_by_rc,
                                  size_t insert_size,
                                  bool use_orientation,
                                  bool handle_Ns,
                                  FileReadFlags flags,
                                  ThreadPool::ThreadPool *pool) {
    return MultifileWrap<TellSeqRead>(tellseq_easy_readers(lib, followed_by_rc, insert_size,
                                                           use_orientation, handle_Ns, flags, pool));
}

ReadStreamList<PairedRead> paired_easy_readers(const SequencingLibraryBase &lib,
                                               bool followed_by_rc,
                                               size_t insert_size,
                                               bool use_orientation,
                                               bool handle_Ns,
                                               FileReadFlags flags,
                                               ThreadPool::ThreadPool *pool) {
    ReadStreamList<PairedRead> streams;
    for (const auto &read_pair : lib.paired_reads()) {
        streams.push_back(PairedEasyStream(read_pair.first, read_pair.second, followed_by_rc, insert_size,
                                           use_orientation, handle_Ns, lib.orientation(), flags, pool));
    }

    for (const auto &read_pair : lib.interlaced_reads()) {
        streams.push_back(PairedEasyStream(read_pair, followed_by_rc, insert_size,
                                           use_orientation, handle_Ns, lib.orientation(), flags, pool));
    }
    return streams;
}

ReadStreamList<TellSeqRead> tellseq_easy_readers(const SequencingLibraryBase &lib,
                                                 bool followed_by_rc,
                                                 size_t insert_size,
                                                 bool use_orientation,
                                                 bool handle_Ns,
                                                 FileReadFlags flags,
                                                 ThreadPool::ThreadPool *pool) {
    VERIFY_MSG(lib.has_aux(), "index stream must be specified for tellseq reads");
    ReadStreamList<TellSeqRead> streams;
    auto aux_it = lib.aux_begin();
    for (const auto &read_pair : lib.paired_reads()) {
        VERIFY(aux_it != lib.aux_end());
        streams.push_back(TellSeqEasyStream(read_pair.first, read_pair.second,
                                            *aux_it++,
                                            followed_by_rc, insert_size,
                                            use_orientation, handle_Ns, lib.orientation(), flags, pool));
    }

    return streams;
}

ReadStreamList<SingleRead> single_easy_readers(const SequencingLibraryBase &lib,
                                               bool followed_by_rc,
                                               bool including_paired_reads,
                                               bool handle_Ns,
                                               FileReadFlags flags,
                                               ThreadPool::ThreadPool *pool) {
    ReadStreamList<SingleRead> streams;
    if (including_paired_reads) {
      for (const auto &read : lib.reads()) {
        //do we need input_file function here?
        streams.push_back(EasyStream(read, followed_by_rc, handle_Ns, flags, pool));
      }
    } else {
      for (const auto &read : lib.single_reads()) {
        streams.push_back(EasyStream(read, followed_by_rc, handle_Ns, flags, pool));
      }
    }
    return streams;
}

SingleStream single_easy_reader(const SequencingLibraryBase &lib,
                                bool followed_by_rc,
                                bool including_paired_reads,
                                bool handle_Ns,
                                FileReadFlags flags,
                                ThreadPool::ThreadPool *pool) {
    return MultifileWrap<io::SingleRead>(
        single_easy_readers(lib, followed_by_rc, including_paired_reads, handle_Ns, flags, pool));
}

ReadStreamList<SingleRead> merged_easy_readers(const SequencingLibraryBase &lib,
                                               bool followed_by_rc,
                                               bool handle_Ns,
                                               FileReadFlags flags,
                                               ThreadPool::ThreadPool *pool) {
    ReadStreamList<SingleRead> streams;
    for (const auto& read : lib.merged_reads()) {
        streams.push_back(EasyStream(read, followed_by_rc, handle_Ns, flags, pool));
    }
    return streams;
}

SingleStream merged_easy_reader(const SequencingLibraryBase &lib,
                                bool followed_by_rc,
                                bool handle_Ns,
                                FileReadFlags flags,
                                ThreadPool::ThreadPool *pool) {
    return MultifileWrap<io::SingleRead>(
        merged_easy_readers(lib, followed_by_rc, handle_Ns, flags, pool));
}

}
