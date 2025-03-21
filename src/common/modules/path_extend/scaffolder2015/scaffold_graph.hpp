
//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

//
// Created by andrey on 17.09.15.
//
#pragma once

#include "utils/logger/logger.hpp"
#include "assembly_graph/core/graph.hpp"
#include "modules/path_extend/paired_library.hpp"
#include "connection_condition2015.hpp"
#include "adt/iterator_range.hpp"

namespace path_extend {
namespace scaffold_graph {

//do NOT add "using namespace debruijn_graph" in order not to confuse between EdgeId typdefs

class ScaffoldGraph {

public:
    //EdgeId in de Bruijn graph is vertex in scaffolding graph
    typedef debruijn_graph::EdgeId ScaffoldVertex;

    //Unique edge id
    typedef size_t ScaffoldEdgeIdT;

    //Scaffold edge indormation class
    struct ScaffoldEdge {
    private:
        //unique id
        ScaffoldEdgeIdT id_;
        //id counter
        static std::atomic<ScaffoldEdgeIdT> scaffold_edge_id_;

        ScaffoldVertex start_;
        ScaffoldVertex end_;
        //color = lib#
        size_t color_;
        //read pair weight or anything else
        double weight_;
        //todo discuss
        size_t length_;

    public:

        ScaffoldEdge(ScaffoldVertex start, ScaffoldVertex end, size_t lib_id = (size_t) -1, double weight = 0, size_t length = 0) :
            id_(scaffold_edge_id_++),
            start_(start), end_(end),
            color_(lib_id),
            weight_(weight),
            length_(length){
        }

        //for consistency with dijkstra
        explicit ScaffoldEdge(size_t ): id_(scaffold_edge_id_++), start_(nullptr), end_(nullptr),
                                        color_((size_t) -1), weight_(0), length_(0) {}

        ScaffoldEdge(): id_(scaffold_edge_id_++), start_(nullptr), end_(nullptr),
                        color_((size_t) -1), weight_(0), length_(0) {}

        ScaffoldEdgeIdT getId() const {
            return id_;
        }


        size_t getColor() const {
            return color_;
        }

        double getWeight() const {
            return weight_;
        }

        size_t getLength() const {
            return length_;
        }

        const ScaffoldVertex getStart() const {
            return start_;
        }

        const ScaffoldVertex getEnd() const {
            return end_;
        }

        bool operator==(const ScaffoldEdge &e) const {
            return color_ == e.color_ && weight_ == e.weight_ && start_ == e.start_ && end_ == e.end_;
        }

        bool operator==(const ScaffoldEdge &e) {
            return color_ == e.color_ && weight_ == e.weight_ && start_ == e.start_ && end_ == e.end_;
        }

        bool operator<(const ScaffoldEdge& rhs) const;
        bool operator>(const ScaffoldEdge& rhs) const;
        bool operator<=(const ScaffoldEdge& rhs) const;
        bool operator>=(const ScaffoldEdge& rhs) const;
    };

    //typedef for possibility to use in templated graph visualizers
    typedef ScaffoldVertex VertexId;
    typedef ScaffoldEdge EdgeId;

    //All vertices are stored in set
    typedef std::set<ScaffoldVertex> VertexStorage;
    //Edges are stored in map: Id -> Edge Information
    typedef std::unordered_map<ScaffoldEdgeIdT, ScaffoldEdge> EdgeStorage;
    //Adjacency list contains vertrx and edge id (instead of whole edge information)
    typedef std::unordered_multimap<ScaffoldVertex, ScaffoldEdgeIdT> AdjacencyStorage;

    struct ConstScaffoldEdgeIterator: public boost::iterator_facade<ConstScaffoldEdgeIterator,
                                                                    const ScaffoldEdge,
                                                                    boost::forward_traversal_tag> {
    private:
        EdgeStorage::const_iterator iter_;

    public:
        ConstScaffoldEdgeIterator(EdgeStorage::const_iterator iter) : iter_(iter) {
        }

    private:
        friend class boost::iterator_core_access;

        void increment() {
            ++iter_;
        }

        bool equal(const ConstScaffoldEdgeIterator &other) const {
            return iter_ == other.iter_;
        }

        const ScaffoldEdge& dereference() const {
            return iter_->second;
        }
    };

//TODO:: fix this. Seems that only ebegin and eend are broken.
private:
    EdgeStorage edges_;

    VertexStorage vertices_;

    const debruijn_graph::Graph &assembly_graph_;

    AdjacencyStorage outgoing_edges_;

    AdjacencyStorage incoming_edges_;

    void AddEdgeSimple(const ScaffoldEdge &e);

    //Delete outgoing edge from adjancecy list without checks
    void DeleteOutgoing(const ScaffoldEdge &e);

    //Delete incoming edge from adjancecy list without checks
    void DeleteIncoming(const ScaffoldEdge &e);

    //Delete all edge info from storage
    void DeleteEdgeFromStorage(const ScaffoldEdge &e);

    //Detelte all outgoing from v edges from  adjacency lists
    void DeleteAllOutgoingEdgesSimple(ScaffoldVertex v);

    //Detelte all incoming from v edges from  adjacency lists
    void DeleteAllIncomingEdgesSimple(ScaffoldVertex v);

public:
    ScaffoldGraph(const debruijn_graph::Graph &g) : assembly_graph_(g) {
    }

    bool Exists(ScaffoldVertex assembly_graph_edge) const;

    bool Exists(const ScaffoldEdge &e) const;

    ScaffoldVertex conjugate(ScaffoldVertex assembly_graph_edge) const;

    //Return structure thay is equal to conjugate of e (not exactrly the same structure as in graph)
    ScaffoldEdge conjugate(const ScaffoldEdge &e) const;

    //Add isolated vertex to the graph if not exitsts
    bool AddVertex(ScaffoldVertex assembly_graph_edge);

    void AddVertices(const std::set<ScaffoldVertex> &vertices);

    //Add edge (and conjugate) if not exists
    //v1 and v2 must exist
    bool AddEdge(ScaffoldVertex v1, ScaffoldVertex v2, size_t lib_id, double weight, size_t length);

    bool AddEdge(const ScaffoldEdge &e);

    //Rempve edge from edge container and all adjacency lists
    bool RemoveEdge(const ScaffoldEdge &e);

    //Remove vertex and all adjacent edges
    bool RemoveVertex(ScaffoldVertex assembly_graph_edge);

    bool IsVertexIsolated(ScaffoldVertex assembly_graph_edge) const;

    VertexStorage::const_iterator vbegin() const;

    VertexStorage::const_iterator vend() const;

    adt::iterator_range<VertexStorage::const_iterator> vertices() const;

    ConstScaffoldEdgeIterator ebegin() const;

    ConstScaffoldEdgeIterator eend() const;

    adt::iterator_range<ScaffoldGraph::ConstScaffoldEdgeIterator> edges() const;

    size_t int_id(ScaffoldVertex v) const;

    size_t int_id(ScaffoldEdge e) const;

    ScaffoldVertex EdgeStart(ScaffoldEdge e) const;

    ScaffoldVertex EdgeEnd(ScaffoldEdge e) const;

    size_t VertexCount() const;

    size_t EdgeCount() const;

    const debruijn_graph::Graph & AssemblyGraph() const;

    std::vector<ScaffoldEdge> OutgoingEdges(ScaffoldVertex assembly_graph_edge) const;

    std::vector<ScaffoldEdge> IncomingEdges(ScaffoldVertex assembly_graph_edge) const;

    size_t OutgoingEdgeCount(ScaffoldVertex assembly_graph_edge) const;

    size_t IncomingEdgeCount(ScaffoldVertex assembly_graph_edge) const;

    bool HasUniqueOutgoing(ScaffoldVertex assembly_graph_edge) const;

    bool HasUniqueIncoming(ScaffoldVertex assembly_graph_edge) const;

    ScaffoldEdge UniqueOutgoing(ScaffoldVertex assembly_graph_edge) const;

    ScaffoldEdge UniqueIncoming(ScaffoldVertex assembly_graph_edge) const;

    void Print(std::ostream &os) const;

    string str(const ScaffoldVertex& vertex) const;

    string str(const ScaffoldGraph::ScaffoldEdge& edge) const;

};

} //scaffold_graph
} //path_extend

