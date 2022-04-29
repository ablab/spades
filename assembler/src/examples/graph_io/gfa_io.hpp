//***************************************************************************
//* Copyright (c) 2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "common/assembly_graph/core/graph.hpp"
#include "common/io/graph/gfa_reader.hpp"
#include "common/io/graph/gfa_writer.hpp"

#include <filesystem>
//TODO: add logging

namespace spades_example {

void SaveToGFA(const debruijn_graph::Graph& g, const std::filesystem::path& save_to);

void ReadFromGFA(debruijn_graph::Graph& g, const std::filesystem::path& read_from);

}
