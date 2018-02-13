#pragma once
#include <common/modules/path_extend/scaffolder2015/scaffold_graph.hpp>
#include "common/adt/iterator_range.hpp"
#include "common/assembly_graph/core/graph.hpp"

namespace contracted_graph {
    class AdjacencyMap {
     public:
        typedef debruijn_graph::VertexId VertexId;
        typedef debruijn_graph::EdgeId EdgeId;
        typedef path_extend::scaffold_graph::ScaffoldGraph::ScaffoldGraphVertex ScaffoldVertex;

     private:
        std::map<debruijn_graph::VertexId, vector<ScaffoldVertex>> data_;

    public:
        typedef std::map<VertexId, vector<ScaffoldVertex>>::const_iterator const_iterator;
        typedef std::map<VertexId, vector<ScaffoldVertex>>::value_type value_type;
        AdjacencyMap() = default;
        AdjacencyMap(const VertexId& vertex, const ScaffoldVertex& edge) : data_({{vertex, {edge}}}) {}
        void InsertPair(const VertexId& vertex, const ScaffoldVertex& edge);

        const_iterator begin() const;
        const_iterator end() const;

    };

    class ContractedGraph {
     public:
        typedef debruijn_graph::VertexId VertexId;
        typedef debruijn_graph::EdgeId EdgeId;
        typedef std::map<VertexId, AdjacencyMap>::const_iterator const_iterator;
        typedef std::set<VertexId>::const_iterator vertex_iterator;
        typedef path_extend::scaffold_graph::ScaffoldGraph::ScaffoldGraphVertex ScaffoldVertex;

     protected:
        std::map<VertexId, AdjacencyMap> outcoming_;
        std::map<VertexId, AdjacencyMap> incoming_;
        std::set<VertexId> vertices_;
        std::map<VertexId, size_t> capacity_;
    public:

        ContractedGraph() = default;
        virtual ~ContractedGraph() = default;

        void InsertVertex(const VertexId& vertex);
        void InsertEdge(const VertexId& head, const VertexId& tail, const ScaffoldVertex& edge);

        AdjacencyMap::const_iterator in_begin(const VertexId& vertex) const;
        AdjacencyMap::const_iterator in_end(const VertexId& vertex) const;

        AdjacencyMap::const_iterator out_begin(const VertexId& vertex) const;
        AdjacencyMap::const_iterator out_end(const VertexId& vertex) const;

        adt::iterator_range<AdjacencyMap::const_iterator> outcoming(const VertexId& vertex) const;
        adt::iterator_range<AdjacencyMap::const_iterator> incoming(const VertexId& vertex) const;

        size_t getOutDegree(const VertexId& vertex) const;
        size_t getInDegree(const VertexId& vertex) const;

        vector <ScaffoldVertex> getIncomingEdges(const VertexId& vertex) const;
        vector <ScaffoldVertex> getOutcomingEdges(const VertexId& vertex) const;

        size_t capacity(const VertexId& vertex) const;
        void InsertCapacity(const VertexId& vertex, size_t capacity);

        bool ContainsVertex(const VertexId& vertex) const;

        vertex_iterator begin() const;
        vertex_iterator end() const;

        size_t size() const;
        size_t CountEdges() const;
    };
}