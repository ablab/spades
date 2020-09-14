//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/graph_support/scaff_supplementary.hpp"
#include "auxiliary_graphs/scaffold_graph/scaffold_graph.hpp"

#include <fstream>

namespace path_extend {
namespace read_cloud {

class ScaffoldGraphStorage {
  public:
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef debruijn_graph::EdgeId EdgeId;

    ScaffoldGraphStorage(ScaffoldGraph &&large_scaffold_graph,
                         ScaffoldGraph &&small_scaffold_graph,
                         const ScaffoldingUniqueEdgeStorage &large_unique_storage,
                         const ScaffoldingUniqueEdgeStorage &small_unique_storage);

    ScaffoldGraphStorage &operator=(const ScaffoldGraphStorage &other);

    const ScaffoldGraph &GetLargeScaffoldGraph() const;
    const ScaffoldGraph &GetSmallScaffoldGraph() const;
    size_t GetLargeLengthThreshold() const;
    size_t GetSmallLengthThreshold() const;
    const ScaffoldingUniqueEdgeStorage &GetLargeUniqueStorage() const;
    const ScaffoldingUniqueEdgeStorage &GetSmallUniqueStorage() const;

    void SetLargeScaffoldGraph(const ScaffoldGraph &large_scaffold_graph);
    void SetSmallScaffoldGraph(const ScaffoldGraph &small_scaffold_graph);

    void Save(const std::string &path) const;
    void Load(const std::string &path, const std::map<size_t, debruijn_graph::EdgeId> &edge_map);

  private:
    void ReplaceScaffoldGraph(const ScaffoldGraph &from, ScaffoldGraph &to);

    ScaffoldGraph large_scaffold_graph_;
    ScaffoldGraph small_scaffold_graph_;
    const ScaffoldingUniqueEdgeStorage &large_unique_storage_;
    const ScaffoldingUniqueEdgeStorage &small_unique_storage_;
    size_t large_length_threshold_;
    size_t small_length_threshold_;
};

class ScaffoldGraphSerializer {
  public:
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef debruijn_graph::EdgeId EdgeId;

    void SaveScaffoldGraph(std::ofstream &fout, const ScaffoldGraph &graph) const;

    void LoadScaffoldGraph(std::ifstream &fin, ScaffoldGraph &graph,
                           const std::map<size_t, debruijn_graph::EdgeId> &edge_map) const;

    DECL_LOGGER("ScaffoldGraphSerializer");
};
}
}