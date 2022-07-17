#pragma once

#include "assembly_graph/paths/bidirectional_path_container.hpp"
#include "assembly_graph/paths/bidirectional_path_io/bidirectional_path_output.hpp"
#include "barcode_index/barcode_info_extractor.hpp"
#include "io/graph/gfa_reader.hpp"

#include "vertex_resolver.hpp"

namespace cont_index {

class PathExtractor {
    typedef std::vector<debruijn_graph::EdgeId> SimplePath;
  public:
    PathExtractor(const debruijn_graph::Graph &graph) : graph_(graph) {}

    void ExtractPaths(path_extend::PathContainer &paths, const VertexResults &vertex_results, bool canonical = true) const;
//    void PrintPaths(const path_extend::PathContainer &paths,
//                    const std::filesystem::path &output_path,
//                    io::IdMapper<std::string> *id_mapper) const;
  private:
    bool IsConjugatePair(const SimplePath &first, const SimplePath &second) const;

    const debruijn_graph::Graph &graph_;
};

}