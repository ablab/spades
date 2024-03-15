// Copyright (c) 2020 Robert Vaser

#ifndef SPOA_GRAPH_HPP_
#define SPOA_GRAPH_HPP_

#include <atomic>
#include <cstdint>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#ifdef SPOA_USE_CEREAL
#include "cereal/access.hpp"
#include "cereal/types/memory.hpp"
#include "cereal/types/vector.hpp"
#include "cereal/types/utility.hpp"
#endif

namespace spoa {

using Alignment = std::vector<std::pair<std::int32_t, std::int32_t>>;

class Graph {
 public:
  Graph();

  Graph(const Graph&) = delete;
  Graph& operator=(const Graph&) = delete;

  Graph(Graph&&) = default;
  Graph& operator=(Graph&&) = default;

  ~Graph() = default;

  struct Node;
  struct Edge;

  struct Node {
   public:
    Node(std::uint32_t id, std::uint32_t code);

    Node(const Node&) = default;
    Node& operator=(const Node&) = default;

    Node(Node&&) = default;
    Node& operator=(Node&&) = default;

    Node* Successor(std::uint32_t label) const;

    std::uint32_t Coverage() const;

    std::uint32_t id;
    std::uint32_t code;
    std::vector<Edge*> inedges;
    std::vector<Edge*> outedges;
    std::vector<Node*> aligned_nodes;

   private:
#ifdef SPOA_USE_CEREAL
    Node() = default;

    template<class Archive>
    void serialize(Archive& archive) {  // NOLINT
      archive(CEREAL_NVP(id), CEREAL_NVP(code));
    }

    friend cereal::access;
#endif
  };
  struct Edge {
   public:
    Edge(Node* tail, Node* head, std::uint32_t label, std::uint32_t weight);

    Edge(const Edge&) = default;
    Edge& operator=(const Edge&) = default;

    Edge(Edge&&) = default;
    Edge& operator=(Edge&&) = default;

    void AddSequence(std::uint32_t label, std::uint32_t weight = 1);

    Node* tail;
    Node* head;
    std::vector<std::uint32_t> labels;
    std::int64_t weight;

   private:
#ifdef SPOA_USE_CEREAL
    Edge() = default;

    template<class Archive>
    void serialize(Archive& archive) {  // NOLINT
      archive(labels, weight);
    }

    friend cereal::access;
#endif
  };

  const std::vector<std::unique_ptr<Node>>& nodes() const {
    return nodes_;
  }

  const std::vector<std::unique_ptr<Edge>>& edges() const {
    return edges_;
  }

  const std::vector<Node*>& rank_to_node() const {
    return rank_to_node_;
  }

  const std::vector<Node*>& sequences() const {
    return sequences_;
  }

  std::uint32_t num_codes() const {
    return num_codes_;
  }

  std::uint8_t coder(std::uint8_t c) const {
    return coder_[c];
  }

  std::uint8_t decoder(std::uint8_t code) const {
    return decoder_[code];
  }

  const std::vector<Node*>& consensus() const {
    return consensus_;
  }

  void AddAlignment(
      const Alignment& alignment,
      const std::string& sequence,
      std::uint32_t weight = 1);

  void AddAlignment(
      const Alignment& alignment,
      const std::string& sequence,
      const std::vector<std::uint32_t>& weights);

  void AddAlignment(
      const Alignment& alignment,
      const std::string& sequence,
      const std::string& quality);

  void AddAlignment(
      const Alignment& alignment,
      const char* sequence, std::uint32_t sequence_len,
      std::uint32_t weight = 1);

  void AddAlignment(
      const Alignment& alignment,
      const char* sequence, std::uint32_t sequence_len,
      const std::vector<std::uint32_t>& weights);

  void AddAlignment(
      const Alignment& alignment,
      const char* sequence, std::uint32_t sequence_len,
      const char* quality, std::uint32_t quality_len);

  std::vector<std::string> GenerateMultipleSequenceAlignment(
      bool include_consensus = false);

  std::string GenerateConsensus();

  std::string GenerateConsensus(
      std::vector<std::uint32_t>* summary,
      bool verbose = false);

  Graph Subgraph(
      std::uint32_t begin,
      std::uint32_t end,
      std::vector<const Node*>* subgraph_to_graph) const;

  void UpdateAlignment(
      const std::vector<const Node*>& subgraph_to_graph,
      Alignment* alignment) const;

  // print with Graphviz
  void PrintDot(const std::string& path) const;

  void Clear();

#ifdef SPOA_USE_CEREAL
  template<class Archive>
  void save(Archive& archive) const {  // NOLINT
    std::vector<std::uint32_t> sequences;
    for (const auto& it : sequences_) {
      sequences.emplace_back(it->id);
    }

    std::vector<std::pair<std::uint32_t, std::uint32_t>> connections;
    for (const auto& it : edges_) {
      connections.emplace_back(it->tail->id, it->head->id);
    }

    std::vector<std::pair<std::uint32_t, std::uint32_t>> aligned_nodes;
    for (const auto& it : nodes_) {
      for (const auto& jt : it->aligned_nodes) {
        if (it->id < jt->id) {
          aligned_nodes.emplace_back(it->id, jt->id);
        }
      }
    }

    std::vector<std::uint32_t> rank_to_node_id;
    for (const auto& it : rank_to_node_) {
      rank_to_node_id.emplace_back(it->id);
    }

    std::vector<std::uint32_t> consensus;
    for (const auto& it : consensus_) {
      consensus.emplace_back(it->id);
    }

    archive(
        num_codes_,
        coder_,
        decoder_,
        nodes_,
        edges_,
        sequences,
        connections,
        aligned_nodes,
        rank_to_node_id,
        consensus);
  }

  template<class Archive>
  void load(Archive& archive) {  // NOLINT
    std::vector<std::uint32_t> sequences;
    std::vector<std::pair<std::uint32_t, std::uint32_t>> connections;
    std::vector<std::pair<std::uint32_t, std::uint32_t>> aligned_nodes;
    std::vector<std::uint32_t> rank_to_node_id;
    std::vector<std::uint32_t> consensus;

    archive(
        num_codes_,
        coder_,
        decoder_,
        nodes_,
        edges_,
        sequences,
        connections,
        aligned_nodes,
        rank_to_node_id,
        consensus);

    for (const auto& it : sequences) {
      sequences_.emplace_back(nodes_[it].get());
    }

    for (std::uint32_t i = 0; i < connections.size(); ++i) {
      edges_[i]->tail = nodes_[connections[i].first].get();
      edges_[i]->head = nodes_[connections[i].second].get();

      edges_[i]->tail->outedges.emplace_back(edges_[i].get());
      edges_[i]->head->inedges.emplace_back(edges_[i].get());
    }

    for (const auto& it : aligned_nodes) {
      nodes_[it.first]->aligned_nodes.emplace_back(nodes_[it.second].get());
      nodes_[it.second]->aligned_nodes.emplace_back(nodes_[it.first].get());
    }

    for (const auto& it : rank_to_node_id) {
      rank_to_node_.emplace_back(nodes_[it].get());
    }

    for (const auto& it : consensus) {
      consensus_.emplace_back(nodes_[it].get());
    }
  }
#endif

 private:
  Node* AddNode(std::uint32_t code);

  void AddEdge(Node* tail, Node* head, std::uint32_t weight);

  Node* AddSequence(
      const char* sequence,
      const std::vector<std::uint32_t>& weights,
      std::uint32_t begin,
      std::uint32_t end);

  void TopologicalSort();

  bool IsTopologicallySorted() const;

  void TraverseHeaviestBundle();

  Node* BranchCompletion(
      std::uint32_t rank,
      std::vector<std::int64_t>* scores,
      std::vector<Node*>* predecessors);

  std::vector<bool> ExtractSubgraph(const Node* begin, const Node* end) const;

  std::vector<std::uint32_t> InitializeMultipleSequenceAlignment(
      std::uint32_t* row_size = nullptr) const;

  std::uint32_t num_codes_;
  std::vector<std::int32_t> coder_;
  std::vector<std::int32_t> decoder_;
  std::vector<Node*> sequences_;
  std::vector<std::unique_ptr<Node>> nodes_;
  std::vector<std::unique_ptr<Edge>> edges_;
  std::vector<Node*> rank_to_node_;
  std::vector<Node*> consensus_;
};

}  // namespace spoa

#endif  // SPOA_GRAPH_HPP_
