//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "logger/logger.hpp"
#include "simple_tools.hpp"
#include "io/io_helper.hpp"
#include "io/library.hpp"

#include "config_struct.hpp"

namespace debruijn_graph {

inline
io::PairedStreamPtr paired_easy_reader(const io::SequencingLibrary<debruijn_config::DataSetData> &lib,
                                       bool followed_by_rc,
                                       size_t insert_size,
                                       bool change_read_order = false,
                                       bool use_orientation = true,
                                       io::OffsetType offset_type = io::PhredOffset) {
  io::ReadStreamList<io::PairedRead> streams;
  for (auto it = lib.paired_begin(); it != lib.paired_end(); ++it) {
      streams.push_back(io::PairedEasyStream(it->first, it->second, followed_by_rc, insert_size, change_read_order,
                                             use_orientation, lib.orientation(), offset_type));
  }
  return io::MultifileWrap<io::PairedRead>(streams);
}

inline
io::SingleStreamPtr single_easy_reader(const io::SequencingLibrary<debruijn_config::DataSetData> &lib,
                                       bool followed_by_rc,
                                       bool including_paired_reads,
                                       io::OffsetType offset_type = io::PhredOffset) {
  io::ReadStreamList<io::SingleRead> streams;
  if (including_paired_reads) {
    for (auto it = lib.reads_begin(); it != lib.reads_end(); ++it) {
      //do we need input_file function here?
      streams.push_back(io::EasyStream(*it, followed_by_rc, offset_type));
    }
  } else {
    for (auto it = lib.single_begin(); it != lib.single_end(); ++it) {
      streams.push_back(io::EasyStream(*it, followed_by_rc, offset_type));
    }
  }

  return io::MultifileWrap<io::SingleRead>(streams);
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
