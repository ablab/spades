//***************************************************************************
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "graph_pack.hpp"

namespace graph_pack {

void FillQuality(GraphPack& gp);
void ClearQuality(GraphPack& gp);

void EnsureIndex(GraphPack& gp);
void EnsureBasicMapping(GraphPack& gp);
void EnsureQuality(GraphPack& gp);
void EnsurePos(GraphPack& gp);
void EnsureDebugInfo(GraphPack& gp);

void InitRRIndices(GraphPack& gp);
void ClearRRIndices(GraphPack& gp);
void ClearPaths(GraphPack& gp);

void DetachAll(GraphPack& gp);
void DetachEdgeIndex(GraphPack& gp);

void PrepareForStage(GraphPack& gp, const char*);

} // namespace debruijn_graph
