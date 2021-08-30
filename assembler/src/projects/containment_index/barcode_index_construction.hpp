//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/core/graph.hpp"
#include "barcode_index/barcode_info_extractor.hpp"
#include "paired_info/paired_info_utils.hpp"
#include "modules/alignment/bwa_sequence_mapper.hpp"

#include "io/binary/read_cloud.hpp"
#include "io/dataset_support/read_converter.hpp"

#include <threadpool/threadpool.hpp>
#include <string>

namespace cont_index {

typedef io::DataSet<debruijn_graph::config::LibraryData> DataSet;
typedef io::SequencingLibrary<debruijn_graph::config::LibraryData> SequencingLib;

void ConstructBarcodeIndex(barcode_index::FrameBarcodeIndex<debruijn_graph::Graph> &barcode_index,
                           paired_info::SequencingLib &lib,
                           const debruijn_graph::Graph &graph,
                           const std::string &workdir,
                           unsigned nthreads,
                           bool bin_load, bool bin_save);

}
