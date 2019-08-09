//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "common/modules/path_extend/scaffolder2015/scaffold_vertex.hpp"
#include "common/adt/iterator_range.hpp"
#include "common/assembly_graph/core/graph.hpp"

namespace contracted_graph {
    class AdjacencyMap {
     public:
        typedef debruijn_graph::VertexId VertexId;
        typedef debruijn_graph::EdgeId EdgeId;
        typedef path_extend::scaffold_graph::ScaffoldVertex ScaffoldVertex;
        typedef std::map<VertexId, std::unordered_set<ScaffoldVertex>>::const_iterator const_iterator;
        typedef std::map<VertexId, std::unordered_set<ScaffoldVertex>>::value_type value_type;

        AdjacencyMap() = default;
        AdjacencyMap(const VertexId& vertex, const ScaffoldVertex& edge) : data_({{vertex, {edge}}}) {}
        void InsertPair(const VertexId& vertex, const ScaffoldVertex& edge);
        void RemovePair(const VertexId &vertex, const ScaffoldVertex &edge);
        bool Contains(const VertexId &vertex, const ScaffoldVertex &edge);

        const_iterator begin() const;
        const_iterator end() const;

      private:
        std::map<debruijn_graph::VertexId, std::unordered_set<ScaffoldVertex>> data_;
    };

    class ContractedGraph {
     public:
        typedef debruijn_graph::VertexId VertexId;
        typedef std::map<VertexId, AdjacencyMap>::const_iterator const_iterator;
        typedef std::set<VertexId>::const_iterator vertex_iterator;
        typedef path_extend::scaffold_graph::ScaffoldVertex ScaffoldVertex;
        typedef ScaffoldVertex EdgeId;

        ContractedGraph(const Graph &assembly_graph_);
        virtual ~ContractedGraph() = default;

        void InsertVertex(const VertexId& vertex);
        void InsertEdge(const VertexId& head, const VertexId& tail, const ScaffoldVertex& edge);
        void RemoveEdge(const VertexId& head, const VertexId& tail, const ScaffoldVertex& edge);
        size_t GetOutDegree(const VertexId &vertex) const;
        size_t GetInDegree(const VertexId &vertex) const;
        std::vector <ScaffoldVertex> GetIncomingEdges(const VertexId &vertex) const;
        std::vector <ScaffoldVertex> GetOutcomingEdges(const VertexId &vertex) const;
        size_t GetCapacity(const VertexId &vertex) const;
        void InsertCapacity(const VertexId& vertex, size_t capacity);
        bool ContainsVertex(const VertexId& vertex) const;

        AdjacencyMap::const_iterator in_begin(const VertexId& vertex) const;
        AdjacencyMap::const_iterator in_end(const VertexId& vertex) const;
        AdjacencyMap::const_iterator out_begin(const VertexId& vertex) const;
        AdjacencyMap::const_iterator out_end(const VertexId& vertex) const;
        adt::iterator_range<AdjacencyMap::const_iterator> outcoming(const VertexId& vertex) const;
        adt::iterator_range<AdjacencyMap::const_iterator> incoming(const VertexId& vertex) const;
        vertex_iterator begin() const;
        vertex_iterator end() const;
        size_t size() const;
        size_t CountEdges() const;

        const Graph &GetAssemblyGraph() const;
        std::string EdgeNucls(EdgeId edge) const;
        std::string VertexNucls(VertexId vertex) const;
        double coverage(EdgeId edge) const;
        size_t length(EdgeId edge) const;
        size_t int_id(EdgeId edge) const;

        ScaffoldVertex conjugate(ScaffoldVertex edge) const;

     protected:
        std::map<VertexId, AdjacencyMap> outcoming_;
        std::map<VertexId, AdjacencyMap> incoming_;
        std::set<VertexId> vertices_;
        std::map<VertexId, size_t> capacity_;

        //for compatibility with visualizers and stuff
        const Graph& assembly_graph_;
    };
}