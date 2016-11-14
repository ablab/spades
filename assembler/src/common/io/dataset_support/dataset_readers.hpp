//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "utils/logger/logger.hpp"
#include "utils/simple_tools.hpp"
#include "io/reads/io_helper.hpp"
#include "pipeline/library.hpp"
#include "pipeline/config_struct.hpp"

namespace io {

inline
PairedStreamPtr paired_easy_reader(const SequencingLibrary<debruijn_graph::config::DataSetData> &lib,
                                   bool followed_by_rc,
                                   size_t insert_size,
                                   bool change_read_order = false,
                                   bool use_orientation = true,
                                   OffsetType offset_type = PhredOffset) {
    ReadStreamList<PairedRead> streams;
    for (auto read_pair : lib.paired_reads()) {
        streams.push_back(PairedEasyStream(read_pair.first, read_pair.second, followed_by_rc, insert_size, change_read_order,
                                           use_orientation, lib.orientation(), offset_type));
    }
    return MultifileWrap<PairedRead>(streams);
}

inline
ReadStreamList<SingleRead> single_easy_readers(const SequencingLibrary<debruijn_graph::config::DataSetData> &lib,
                                               bool followed_by_rc,
                                               bool including_paired_reads,
                                               bool handle_Ns = true,
                                               OffsetType offset_type = PhredOffset) {
    ReadStreamList<SingleRead> streams;
    if (including_paired_reads) {
      for (const auto& read : lib.reads()) {
        //do we need input_file function here?
        streams.push_back(EasyStream(read, followed_by_rc, handle_Ns, offset_type));
      }
    } else {
      for (const auto& read : lib.single_reads()) {
        streams.push_back(EasyStream(read, followed_by_rc, handle_Ns, offset_type));
      }
    }
    return streams;
}

inline
SingleStreamPtr single_easy_reader(const SequencingLibrary<debruijn_graph::config::DataSetData> &lib,
                                   bool followed_by_rc,
                                   bool including_paired_reads,
                                   bool handle_Ns = true,
                                   OffsetType offset_type = PhredOffset) {
    return MultifileWrap<io::SingleRead>(
           single_easy_readers(lib, followed_by_rc, including_paired_reads, handle_Ns, offset_type));
}

inline
PairedStreamPtr paired_easy_reader_for_libs(std::vector<size_t> libs,
                                            bool followed_by_rc,
                                            size_t insert_size,
                                            bool change_read_order = false,
                                            bool use_orientation = true,
                                            OffsetType offset_type = PhredOffset) {
    ReadStreamList<io::PairedRead> streams;
    for (size_t i = 0; i < libs.size(); ++i) {
      streams.push_back(paired_easy_reader(cfg::get().ds.reads[libs[i]],
                                           followed_by_rc, insert_size, change_read_order, use_orientation, offset_type));
    }
    return MultifileWrap<PairedRead>(streams);
}


inline
PairedStreamPtr paired_easy_reader(bool followed_by_rc,
                                   size_t insert_size,
                                   bool change_read_order = false,
                                   bool use_orientation = true,
                                   OffsetType offset_type = PhredOffset) {

    std::vector<size_t> all_libs(cfg::get().ds.reads.lib_count());
    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i)
        all_libs[i] = i;

    // FIXME: Should we use only first library?
    // No, this one is for all libs together
    return paired_easy_reader_for_libs(all_libs, followed_by_rc, insert_size, change_read_order, use_orientation, offset_type);
}


inline
SingleStreamPtr single_easy_reader_for_libs(vector<size_t> libs,
                                            bool followed_by_rc,
                                            bool including_paired_reads,
                                            OffsetType offset_type = PhredOffset) {
    ReadStreamList<SingleRead> streams;
    for (size_t i = 0; i < libs.size(); ++i) {
        streams.push_back(single_easy_reader(cfg::get().ds.reads[libs[i]],
                                             followed_by_rc, including_paired_reads, offset_type));
    }
    return MultifileWrap<SingleRead>(streams);
}

inline
SingleStreamPtr single_easy_reader(bool followed_by_rc,
                                   bool including_paired_reads,
                                   OffsetType offset_type = PhredOffset) {

    std::vector<size_t> all_libs(cfg::get().ds.reads.lib_count());
    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i)
        all_libs[i] = i;

    return single_easy_reader_for_libs(all_libs, followed_by_rc, including_paired_reads, offset_type);
}

}
