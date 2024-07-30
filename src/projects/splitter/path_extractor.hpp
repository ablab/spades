//***************************************************************************
//* Copyright (c) 2021-2023 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/paths/bidirectional_path_container.hpp"
#include "assembly_graph/paths/bidirectional_path_io/bidirectional_path_output.hpp"
#include "barcode_index/barcode_info_extractor.hpp"
#include "io/graph/gfa_reader.hpp"

#include "vertex_resolver.hpp"

namespace cont_index {

class PathExtractor {
    typedef std::vector<debruijn_graph::EdgeId> SimplePath;
    typedef std::unordered_map<debruijn_graph::EdgeId, std::unordered_set<debruijn_graph::EdgeId>> VertexLinkStorage;
  public:
    explicit PathExtractor(const debruijn_graph::Graph &graph) : graph_(graph) {}

    void ExtractPaths(path_extend::PathContainer &paths, const VertexResults &vertex_results) const;
  private:
    struct ScaffoldLinks {
      typedef std::unordered_map<debruijn_graph::EdgeId, size_t> DegreeMap;
      DegreeMap in_degrees;
      DegreeMap out_degrees;
      std::unordered_map<debruijn_graph::EdgeId, debruijn_graph::EdgeId> in_to_out;
      std::unordered_map<debruijn_graph::EdgeId, std::unordered_set<debruijn_graph::EdgeId>> vertex_link_storage;

      ScaffoldLinks(const DegreeMap &in_degrees,
                    const DegreeMap &out_degrees,
                    const std::unordered_map<debruijn_graph::EdgeId, debruijn_graph::EdgeId> &in_to_out,
                    const std::unordered_map<debruijn_graph::EdgeId,
                                             std::unordered_set<debruijn_graph::EdgeId>> &vertex_link_storage)
          : in_degrees(in_degrees),
            out_degrees(out_degrees),
            in_to_out(in_to_out),
            vertex_link_storage(vertex_link_storage) {}
    };

    bool IsConjugatePair(const SimplePath &first, const SimplePath &second) const;
    bool IsGraphLink(const debruijn_graph::EdgeId &first,
                     const debruijn_graph::EdgeId &second,
                     const VertexLinkStorage &vertex_storage) const;

    ScaffoldLinks GetScaffoldLinks(const VertexResults &vertex_results) const;

    const debruijn_graph::Graph &graph_;

    DECL_LOGGER("PathExtractor");
};

}