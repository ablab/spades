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

namespace debruijn_graph {

inline
io::PairedStreamPtr paired_easy_reader(const io::SequencingLibrary<config::DataSetData> &lib,
                                       bool followed_by_rc,
                                       size_t insert_size,
                                       bool change_read_order = false,
                                       bool use_orientation = true,
                                       io::OffsetType offset_type = io::PhredOffset) {
  io::ReadStreamList<io::PairedRead> streams;
  for (auto read_pair : lib.paired_reads()) {
      streams.push_back(io::PairedEasyStream(read_pair.first, read_pair.second, followed_by_rc, insert_size, change_read_order,
                                             use_orientation, lib.orientation(), offset_type));
  }
  return io::MultifileWrap<io::PairedRead>(streams);
}

inline
io::ReadStreamList<io::SingleRead> single_easy_readers(const io::SequencingLibrary<config::DataSetData> &lib,
                                       bool followed_by_rc,
                                       bool including_paired_reads,
                                       bool handle_Ns = true,
                                       io::OffsetType offset_type = io::PhredOffset) {
  io::ReadStreamList<io::SingleRead> streams;
  if (including_paired_reads) {
    for (const auto& read : lib.reads()) {
      //do we need input_file function here?
      streams.push_back(io::EasyStream(read, followed_by_rc, handle_Ns, offset_type));
    }
  } else {
    for (const auto& read : lib.single_reads()) {
      streams.push_back(io::EasyStream(read, followed_by_rc, handle_Ns, offset_type));
    }
  }
  return streams;
}

inline
io::SingleStreamPtr single_easy_reader(const io::SequencingLibrary<config::DataSetData> &lib,
                                       bool followed_by_rc,
                                       bool including_paired_reads,
                                       bool handle_Ns = true,
                                       io::OffsetType offset_type = io::PhredOffset) {
  return io::MultifileWrap<io::SingleRead>(
          single_easy_readers(lib, followed_by_rc, including_paired_reads, handle_Ns, offset_type));
}

inline
io::PairedStreamPtr paired_easy_reader_for_libs(std::vector<size_t> libs,
                                                bool followed_by_rc,
                                                size_t insert_size,
                                                bool change_read_order = false,
                                                bool use_orientation = true,
                                                io::OffsetType offset_type = io::PhredOffset) {
  io::ReadStreamList<io::PairedRead> streams;
  for (size_t i = 0; i < libs.size(); ++i) {
    streams.push_back(paired_easy_reader(cfg::get().ds.reads[libs[i]],
                                         followed_by_rc, insert_size, change_read_order, use_orientation, offset_type));
  }
  return io::MultifileWrap<io::PairedRead>(streams);
}


inline
io::PairedStreamPtr paired_easy_reader(bool followed_by_rc,
                                       size_t insert_size,
                                       bool change_read_order = false,
                                       bool use_orientation = true,
                                       io::OffsetType offset_type = io::PhredOffset) {

  std::vector<size_t> all_libs(cfg::get().ds.reads.lib_count());
  for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i)
    all_libs[i] = i;

  // FIXME: Should we use only first library?
  // No, this one is for all libs together
  return paired_easy_reader_for_libs(all_libs, followed_by_rc, insert_size, change_read_order, use_orientation, offset_type);
}


inline
io::SingleStreamPtr single_easy_reader_for_libs(vector<size_t> libs,
                                                bool followed_by_rc,
                                                bool including_paired_reads,
                                                io::OffsetType offset_type = io::PhredOffset) {
  io::ReadStreamList<io::SingleRead> streams;
  for (size_t i = 0; i < libs.size(); ++i) {
    streams.push_back(single_easy_reader(cfg::get().ds.reads[libs[i]],
                                         followed_by_rc, including_paired_reads, offset_type));
  }
  return io::MultifileWrap<io::SingleRead>(streams);
}

inline
io::SingleStreamPtr single_easy_reader(bool followed_by_rc,
                                       bool including_paired_reads,
                                       io::OffsetType offset_type = io::PhredOffset) {

  std::vector<size_t> all_libs(cfg::get().ds.reads.lib_count());
  for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i)
    all_libs[i] = i;

  return single_easy_reader_for_libs(all_libs, followed_by_rc, including_paired_reads, offset_type);
}

}
