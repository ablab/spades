//***************************************************************************
//* Copyright (c) 2023-2025 SPAdes team
//* Copyright (c) 2019-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "auxiliary_graphs/scaffold_graph/scaffold_vertex.hpp"
#include "adt/iterator_range.hpp"
#include "assembly_graph/core/graph.hpp"

namespace contracted_graph {
class AdjacencyMap {
 public:
    typedef debruijn_graph::VertexId VertexId;
    typedef debruijn_graph::EdgeId EdgeId;
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef std::map<VertexId, std::unordered_set<ScaffoldVertex>> ContractedMap;
    typedef ContractedMap::const_iterator const_iterator;
    typedef ContractedMap::value_type value_type;
    typedef ContractedMap::mapped_type mapped_type;

    AdjacencyMap() = default;
    AdjacencyMap(VertexId vertex, const ScaffoldVertex &edge) : data_({{vertex, {edge}}}) {}
    void InsertPair(VertexId vertex, const ScaffoldVertex &edge);
    void RemovePair(VertexId vertex, const ScaffoldVertex &edge);
    bool Contains(VertexId vertex, const ScaffoldVertex &edge);
    bool empty() const;
    size_t size() const;

    const_iterator begin() const;
    const_iterator end() const;

  private:
    std::map<debruijn_graph::VertexId, std::unordered_set<ScaffoldVertex>> data_;
};

class ContractedGraph {
 public:
    typedef debruijn_graph::VertexId VertexId;
    typedef debruijn_graph::Graph Graph;
    typedef std::set<VertexId> VertexContainer;
    typedef std::map<VertexId, AdjacencyMap> EdgeContainer;
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef AdjacencyMap::const_iterator const_entry_iterator;
    typedef VertexContainer::const_iterator const_vertex_iterator;
    typedef AdjacencyMap::mapped_type::const_iterator internal_edge_iterator;
    typedef ScaffoldVertex EdgeId;

    class const_edge_iterator : public boost::iterator_facade<const_edge_iterator,
                                                              const ScaffoldVertex,
                                                              boost::forward_traversal_tag> {
      public:
        explicit const_edge_iterator(const_entry_iterator entry_it,
                                     internal_edge_iterator edge_it,
                                     const_entry_iterator entry_end):
            entry_it_(entry_it),
            edge_it_(edge_it),
            entry_end_(entry_end) {}

      private:
        friend class boost::iterator_core_access;

        bool equal(const const_edge_iterator &other) const {
            return entry_it_ == other.entry_it_ and edge_it_ == other.edge_it_ and entry_end_ == other.entry_end_;
        }

        const ScaffoldVertex &dereference() const {
            return *edge_it_;
        }

        void increment() {
            ++edge_it_;
            if (edge_it_ == entry_it_->second.end()) {
                ++entry_it_;
                if (entry_it_ != entry_end_) {
                    edge_it_ = entry_it_->second.begin();
                }
            }
        }

        const_entry_iterator entry_it_;
        internal_edge_iterator edge_it_;
        const_entry_iterator entry_end_;
    };

    explicit ContractedGraph(const Graph &assembly_graph);
    virtual ~ContractedGraph() = default;
    ContractedGraph(ContractedGraph &&other) = default;

    void InsertVertex(VertexId vertex);
    void InsertEdge(VertexId head, VertexId tail, const ScaffoldVertex &edge);
    void RemoveEdge(VertexId head, VertexId tail, const ScaffoldVertex &edge);
    size_t GetOutDegree(VertexId vertex) const;
    size_t GetInDegree(VertexId vertex) const;
    size_t GetCapacity(VertexId vertex) const;
    void InsertCapacity(VertexId vertex, size_t capacity);
    bool ContainsVertex(VertexId vertex) const;

    const_entry_iterator in_entry_begin(VertexId vertex) const;
    const_entry_iterator in_entry_end(VertexId vertex) const;
    adt::iterator_range<const_entry_iterator> IncomingEntries(VertexId vertex) const;
    const_entry_iterator out_entry_begin(VertexId vertex) const;
    const_entry_iterator out_entry_end(VertexId vertex) const;
    adt::iterator_range<const_entry_iterator> OutcomingEntries(VertexId vertex) const;

    const_edge_iterator in_edge_begin(VertexId vertex) const;
    const_edge_iterator in_edge_end(VertexId vertex) const;
    adt::iterator_range<const_edge_iterator> IncomingEdges(VertexId vertex) const;
    size_t IncomingEdgeCount(VertexId vertex) const;
    const_edge_iterator out_edge_begin(VertexId vertex) const;
    const_edge_iterator out_edge_end(VertexId vertex) const;
    adt::iterator_range<const_edge_iterator> OutgoingEdges(VertexId vertex) const;
    size_t OutgoingEdgeCount(VertexId vertex) const;

    const_vertex_iterator begin() const;
    const_vertex_iterator end() const;
    adt::iterator_range<const_vertex_iterator> vertices() const;
    size_t size() const;
    size_t CountEdges() const;

    //fixme also iterates over short edges
    auto canonical_edges () const;

    const Graph &GetAssemblyGraph() const;
//    std::string EdgeNucls(EdgeId edge) const;
//    std::string VertexNucls(VertexId vertex) const;
    Sequence EdgeNucls(EdgeId edge) const;
    double coverage(EdgeId edge) const;
    size_t length(EdgeId edge) const;
    size_t int_id(EdgeId edge) const;

    ScaffoldVertex conjugate(ScaffoldVertex edge) const;
    VertexId conjugate(VertexId vertex) const;

 protected:
    EdgeContainer outcoming_;
    EdgeContainer incoming_;
    VertexContainer vertices_;
    std::map<VertexId, size_t> capacity_;

    //for edge iterator
    std::unordered_set<ScaffoldVertex> empty_;

    //for compatibility with visualizers and stuff
    const Graph &assembly_graph_;
    };
}