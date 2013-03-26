//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "logger/logger.hpp"
#include <simple_tools.hpp>
#include <io/easy_reader.hpp>
#include <io/single_read.hpp>

namespace debruijn_graph {

typedef io::IReader<io::SingleRead> ReadStream;
typedef io::IReader<io::PairedRead> PairedReadStream;
typedef io::MultifileReader<io::PairedRead> MultiPairedStream;
typedef io::MultifileReader<io::SingleRead> MultiSingleStream;


std::auto_ptr<PairedReadStream> paired_easy_reader_for_libs(vector<size_t> libs,
                                                   bool followed_by_rc,
                                                   size_t insert_size,
                                                   bool change_read_order = false,
                                                   bool revert_second = true,
                                                   io::OffsetType offset_type = io::PhredOffset) {
  std::vector<PairedReadStream*> streams;
  for (size_t i = 0; i < libs.size(); ++i) {
      std::auto_ptr<PairedReadStream> reader = cfg::get().ds.reads[libs[i]].paired_easy_reader(followed_by_rc, insert_size, change_read_order, revert_second, offset_type);
      streams.push_back(reader.get());
      reader.release();
  }
  return auto_ptr<PairedReadStream>(new MultiPairedStream(streams, true));
}


std::auto_ptr<PairedReadStream> paired_easy_reader(bool followed_by_rc,
                                                   size_t insert_size,
                                                   bool change_read_order = false,
                                                   bool revert_second = true,
                                                   io::OffsetType offset_type = io::PhredOffset) {

  std::vector<size_t> all_libs(cfg::get().ds.reads.lib_count());
  for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i)
      all_libs[i] = i;

  // FIXME: Should we use only first library?
  // No, this one is for all libs together
  return paired_easy_reader_for_libs(all_libs, followed_by_rc, insert_size, change_read_order, revert_second, offset_type);
}


std::auto_ptr<ReadStream> single_easy_reader_for_libs(vector<size_t> libs,
                                                    bool followed_by_rc,
                                                    bool including_paired_reads,
                                                    io::OffsetType offset_type = io::PhredOffset) {
  std::vector<ReadStream*> streams;
  for (size_t i = 0; i < libs.size(); ++i) {
      std::auto_ptr<ReadStream> reader = cfg::get().ds.reads[libs[i]].single_easy_reader(followed_by_rc, including_paired_reads, offset_type);
      streams.push_back(reader.get());
      reader.release();
  }
  return auto_ptr<ReadStream>(new MultiSingleStream(streams, true));
}


std::auto_ptr<ReadStream> single_easy_reader(bool followed_by_rc,
                                             bool including_paired_reads,
                                             io::OffsetType offset_type = io::PhredOffset) {

  std::vector<size_t> all_libs(cfg::get().ds.reads.lib_count());
  for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i)
    all_libs[i] = i;

  return single_easy_reader_for_libs(all_libs, followed_by_rc, including_paired_reads, offset_type);
}

}
