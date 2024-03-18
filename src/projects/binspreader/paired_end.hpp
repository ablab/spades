//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2021-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/core/graph.hpp"

#include "library/library.hpp"
#include "library/library_data.hpp"

#include <filesystem>

typedef io::DataSet<debruijn_graph::config::LibraryData> DataSet;
typedef io::SequencingLibrary<debruijn_graph::config::LibraryData> SequencingLib;

namespace binning {
class LinkIndex;

void FillPairedEndLinks(LinkIndex &pe_links,
                        SequencingLib &lib,
                        const debruijn_graph::Graph &graph,
                        const std::filesystem::path &workdir,
                        unsigned nthreads,
                        bool bin_load, bool bin_save);
};
