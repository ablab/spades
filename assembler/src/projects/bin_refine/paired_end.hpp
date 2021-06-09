//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/core/graph.hpp"

#include "pipeline/library.hpp"
#include "pipeline/library_data.hpp"

typedef io::DataSet<debruijn_graph::config::LibraryData> DataSet;
typedef io::SequencingLibrary<debruijn_graph::config::LibraryData> SequencingLib;

namespace binning {
class LinkIndex;

void FillPairedEndLinks(LinkIndex &pe_links,
                        SequencingLib &lib,
                        const debruijn_graph::Graph &graph,
                        const std::string &workdir,
                        unsigned nthreads,
                        bool bin_load, bool bin_save);
};
