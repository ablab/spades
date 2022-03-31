#pragma once

#include "barcode_index/barcode_info_extractor.hpp"

#include "io/graph/gfa_reader.hpp"
#include "scaffold_graph_helper.hpp"

namespace cont_index {

enum class VertexState {
    Completely,
    Partially,
    Ambiguous,
    Uncovered
};

struct VertexResult {
    VertexResult(VertexState state,
                 const size_t &total_links,
                 const size_t &supporting_links,
                 const std::unordered_map<debruijn_graph::EdgeId, debruijn_graph::EdgeId> &supported_pairs) :
                    state(state),
                    total_links(total_links),
                    supporting_links(supporting_links),
                    supported_pairs(supported_pairs) {}

    VertexState state;
    size_t total_links;
    size_t supporting_links;
    std::unordered_map<debruijn_graph::EdgeId, debruijn_graph::EdgeId> supported_pairs;
};

struct VertexResults {
  VertexResults(const std::unordered_map<debruijn_graph::VertexId, VertexResult> &vertex_to_result) :
      vertex_to_result(vertex_to_result) {}

  std::unordered_map<debruijn_graph::VertexId, VertexResult> vertex_to_result;
};

class PathExtractor {
    typedef std::vector<debruijn_graph::EdgeId> SimplePath;
  public:
    PathExtractor(const debruijn_graph::Graph &graph) : graph_(graph) {}

    std::vector<SimplePath> ExtractPaths(const VertexResults &vertex_results, bool canonical = true) const;
    void PrintPaths(const std::vector<SimplePath> &paths,
                    const std::string &output_path,
                    io::IdMapper<std::string> *id_mapper) const;
  private:
    bool IsConjugatePair(const SimplePath &first, const SimplePath &second) const;

    const debruijn_graph::Graph &graph_;
};

class VertexResolver {
  public:
    typedef std::unordered_map<debruijn_graph::EdgeId, VertexState> ResolutionResults;
    typedef debruijn_graph::EdgeId EdgeId;

    VertexResolver(debruijn_graph::Graph &graph,
                   const std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> &barcode_extractor_ptr,
                   size_t count_threshold,
                   size_t tail_threshold,
                   size_t length_threshold,
                   size_t threads,
                   double score_threshold,
                   const string &output_path);

    VertexResults ResolveVertices();
    VertexResult ResolveVertex(const debruijn_graph::VertexId &vertex,
                               LinkIndexGraphConstructor::BarcodeScoreFunctionPtr score_function) const;

    void PrintVertexResults(const VertexResults &results,
                            const std::string &output_path,
                            io::IdMapper<std::string> *id_mapper) const;
    std::string VertexResultString(const debruijn_graph::VertexId &vertex,
                                   const VertexResult &vertex_result,
                                   io::IdMapper<std::string> *id_mapper) const;
  private:
    VertexState GetState(const std::unordered_map<EdgeId, EdgeId> &in_to_out,
                         const debruijn_graph::VertexId &vertex,
                         bool is_ambiguous,
                         bool is_covered) const;

    debruijn_graph::Graph &graph_;
    std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr_;
    size_t count_threshold_;
    size_t tail_threshold_;
    size_t length_threshold_;
    size_t threads_;
    double score_threshold_;
    const std::string &output_path_;
};

}