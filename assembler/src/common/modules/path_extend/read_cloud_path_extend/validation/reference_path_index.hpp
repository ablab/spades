//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "transition_extractor.hpp"
#include "assembly_graph/core/graph.hpp"

namespace path_extend {
namespace read_cloud {
namespace validation {

class ReferencePathIndex {
    typedef debruijn_graph::EdgeId EdgeId;
    struct EdgeInfo {
      //Position of the edge in the path
      size_t edge_pos_;
      size_t conj_edge_pos_;

      //Position of the start\end vertex in the genome
      size_t start_pos_;
      size_t end_pos_;

      //Genome\chromosome id
      size_t path_;

      EdgeInfo(size_t edge_pos, size_t conj_edge_pos, size_t start_pos, size_t end_pos, size_t path) :
          edge_pos_(edge_pos), conj_edge_pos_(conj_edge_pos), start_pos_(start_pos), end_pos_(end_pos), path_(path) {}
    };

    std::unordered_map<EdgeId, EdgeInfo> edge_to_info_;

  public:
    void Insert(EdgeId edge, size_t edge_pos, size_t conj_edge_pos, size_t start_pos, size_t end_pos, size_t path);
    EdgeInfo at(const EdgeId &edge) const;
    bool Contains(const EdgeId &edge) const;
};

class ReferencePathIndexBuilder {
  public:
    typedef std::vector<std::vector<EdgeWithMapping>> ReferencePaths;

    ReferencePathIndex BuildReferencePathIndex(const UniqueReferencePaths &unique_paths);
};
}
}
}