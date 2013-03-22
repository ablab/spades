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

std::auto_ptr<PairedReadStream> paired_easy_reader(bool followed_by_rc,
                                                   size_t insert_size,
                                                   bool change_read_order = false,
                                                   bool revert_second = true,
                                                   io::OffsetType offset_type = io::PhredOffset) {
  std::vector<PairedReadStream*> streams;
  const io::DataSet &dataset = cfg::get().ds.dataset;
  // FIXME: Should we use only first library?
  for (auto it = dataset.paired_begin(); it != dataset.paired_end(); ++it) {
    io::PairedEasyReader* reader = new io::PairedEasyReader(*it, followed_by_rc, insert_size, change_read_order, revert_second, offset_type);
    streams.push_back(reader);
  }
  return auto_ptr<PairedReadStream>(new MultiPairedStream(streams, true));
}

std::auto_ptr<ReadStream> single_easy_reader(bool followed_by_rc,
                                             bool including_paired_reads,
                                             io::OffsetType offset_type = io::PhredOffset) {
  std::vector<ReadStream*> streams;
  const io::DataSet &dataset = cfg::get().ds.dataset;
  // FIXME: Should we use only first library?
  if (including_paired_reads) {
    for (auto it = dataset.reads_begin(); it != dataset.reads_end(); ++it) {
      streams.push_back(new io::EasyReader(input_file(*it), followed_by_rc, offset_type));
      DEBUG("Using input file: " << input_file(*it));
    }
  } else {
    for (auto it = dataset.single_begin(); it != dataset.single_end(); ++it) {
      streams.push_back(new io::EasyReader(input_file(*it), followed_by_rc, offset_type));
      DEBUG("Using input file: " << input_file(*it));
    }
  }

  return std::auto_ptr<ReadStream>(new MultiSingleStream(streams, true));
}

}
