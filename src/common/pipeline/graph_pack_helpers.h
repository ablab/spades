//***************************************************************************
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "graph_pack.hpp"

namespace graph_pack {

void FillQuality(graph_pack::GraphPack& gp);
void ClearQuality(graph_pack::GraphPack& gp);

void EnsureIndex(graph_pack::GraphPack& gp);
void EnsureBasicMapping(graph_pack::GraphPack& gp);
void EnsureQuality(graph_pack::GraphPack& gp);
void EnsurePos(graph_pack::GraphPack& gp);
void EnsureDebugInfo(graph_pack::GraphPack& gp);

void InitRRIndices(graph_pack::GraphPack& gp);
void ClearRRIndicesAndPaths(graph_pack::GraphPack& gp);

void DetachAll(graph_pack::GraphPack& gp);
void DetachEdgeIndex(graph_pack::GraphPack& gp);

void PrepareForStage(graph_pack::GraphPack& gp, const char*);

} // namespace debruijn_graph
