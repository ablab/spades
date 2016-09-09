//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "utils/simple_tools.hpp"
#include "math/xmath.h"
#include "pipeline/config_struct.hpp"
#include "assembly_graph/paths/mapping_path.hpp"

#include <iostream>
#include <fstream>
#include <map>

namespace debruijn_graph {
namespace stats {

using namespace math;
using namespace omnigraph;

class AbstractStatCounter {
public:
    AbstractStatCounter() {
    }

    virtual ~AbstractStatCounter() {
    }

    virtual void Count() = 0;
    //protected:
    //  DECL_LOGGER("StatCounter")
};

class StatList : AbstractStatCounter {
private:
    vector<AbstractStatCounter *> to_count_;
public:
    StatList(vector<AbstractStatCounter *> to_count =
    vector<AbstractStatCounter *>()) :
            to_count_(to_count) {
    }

    virtual ~StatList() {
    }

    void AddStat(AbstractStatCounter *new_stat) {
        to_count_.push_back(new_stat);
    }

    const vector<AbstractStatCounter *> stats() {
        return to_count_;
    }

    virtual void Count() {
        for (size_t i = 0; i < to_count_.size(); i++) {
            to_count_[i]->Count();
        }
    }

    void DeleteStats() {
        for (size_t i = 0; i < to_count_.size(); i++)
            delete to_count_[i];
        to_count_.clear();
    }
};

template<class Graph>
class VertexEdgeStat : public AbstractStatCounter {
private:
    const Graph &graph_;
public:
    VertexEdgeStat(const Graph &graph) :
            graph_(graph) {
    }

    virtual ~VertexEdgeStat() {
    }

    size_t vertices() {
        return graph_.size();
    }

    size_t edges() {
        size_t edgeNumber = 0;
        size_t sum_edge_length = 0;
        for (auto iterator = graph_.ConstEdgeBegin(); !iterator.IsEnd();
             ++iterator) {
            edgeNumber++;
            //      if (graph_.coverage(*iterator) > 30) {
            sum_edge_length += graph_.length(*iterator);
            //      }
        }
        return edgeNumber;
    }

    size_t edge_length() {
        size_t sum_edge_length = 0;
        for (auto iterator = graph_.ConstEdgeBegin(); !iterator.IsEnd();
             ++iterator) {
            if (graph_.coverage(*iterator) > 30) {
                sum_edge_length += graph_.length(*iterator);
            }
        }
        return sum_edge_length;
    }

    virtual void Count() {
        INFO(
                "Vertex count=" << vertices() << "; Edge count=" << edges());
        INFO(
                "sum length of edges " << edge_length());
    }
};

template<class Graph>
class BlackEdgesStat : public AbstractStatCounter {
private:
    typedef typename Graph::EdgeId EdgeId;
    const Graph &graph_;
    Path<EdgeId> path1_;
    Path<EdgeId> path2_;
public:
    BlackEdgesStat(const Graph &graph, Path<EdgeId> path1, Path<EdgeId> path2) :
            graph_(graph), path1_(path1), path2_(path2) {
    }

    virtual ~BlackEdgesStat() {
    }

    virtual void Count() {
        size_t black_count = 0;
        size_t edge_count = 0;
        const vector <EdgeId> path_edges1 = path1_.sequence();
        const vector <EdgeId> path_edges2 = path2_.sequence();
        set <EdgeId> colored_edges;
        colored_edges.insert(path_edges1.begin(), path_edges1.end());
        colored_edges.insert(path_edges2.begin(), path_edges2.end());
        size_t sum_length = 0;
        for (auto it = graph_.ConstEdgeBegin(); !it.IsEnd(); ++it) {
            edge_count++;
            if (colored_edges.count(*it) == 0) {
                black_count++;
                sum_length += graph_.length(*it);
            }
        }
        if (edge_count > 0) {
            INFO("Error edges count: " << black_count << " which is " <<
                 100.0 * (double) black_count / (double) edge_count << "% of all edges");
            INFO("Total length of all black edges: " << sum_length << ". While double genome length is " <<
                 (2 * cfg::get().ds.reference_genome.size()));
        } else {
            INFO("Error edges count: " << black_count << " which is 0% of all edges");
        }
    }
};

template<class Graph>
class NStat : public AbstractStatCounter {
private:
    typedef typename Graph::EdgeId EdgeId;
    const Graph &graph_;
    Path<EdgeId> path_;
    size_t perc_;
public:
    NStat(const Graph &graph, Path<EdgeId> path, size_t perc = 50) :
            graph_(graph), path_(path), perc_(perc) {
    }

    virtual ~NStat() {
    }

    virtual void Count() {
        vector <size_t> lengths;
        size_t sum_all = 0;
        for (size_t i = 0; i < path_.size(); i++) {
            lengths.push_back(graph_.length(path_[i]));
            sum_all += graph_.length(path_[i]);
        }
        sort(lengths.begin(), lengths.end());
        size_t sum = 0;
        size_t current = lengths.size();
        while (current > 0 && (double) sum < (double) perc_ * 0.01 * (double) sum_all) {
            current--;
            sum += lengths[current];
        }
        if (current < lengths.size())
            INFO("N" << perc_ << ": " << lengths[current]);
    }
};

template<class Graph>
class IsolatedEdgesStat : public AbstractStatCounter {
private:
    typedef typename Graph::EdgeId EdgeId;
    const Graph &graph_;
    set <EdgeId> black_edges_;
    vector <size_t> lengths;
public:
    IsolatedEdgesStat(const Graph &graph, Path<EdgeId> path1,
                      Path<EdgeId> path2) :
            graph_(graph) {
        for (auto it = graph.ConstEdgeBegin(); !it.IsEnd(); ++it) {
            black_edges_.insert(*it);
        }
        for (size_t i = 0; i < path1.size(); i++) {
            black_edges_.erase(path1[i]);
        }
        for (size_t i = 0; i < path2.size(); i++) {
            black_edges_.erase(path2[i]);
        }
    }

    virtual ~IsolatedEdgesStat() {
    }

    virtual void Count() {
        lengths.clear();
        for (auto it = graph_.ConstEdgeBegin(); !it.IsEnd(); ++it) {
            EdgeId edge = *it;
            if (graph_.IsDeadEnd(graph_.EdgeEnd(edge))
                && graph_.IsDeadStart(graph_.EdgeStart(edge))
                && black_edges_.count(edge) == 0) {
                lengths.push_back(graph_.length(edge));
            }
        }
        INFO("Isolated not black edges: " << lengths.size());
        WriteLengths(cfg::get().output_dir, "isolated_edges.txt");
    }

    void WriteLengths(string folder_name, string file_name) {
        ofstream os;
        os.open((folder_name + "/" + file_name).c_str());
        WriteLengths(os);
        os.close();
    }

    void WriteLengths(ostream &os) {
        sort(lengths.begin(), lengths.end());
        for (size_t i = 0; i < lengths.size(); i++) {
            os << lengths[i] << endl;
        }
    }
};

template<class Graph>
class SelfComplementStat : public AbstractStatCounter {
private:
    typedef typename Graph::EdgeId EdgeId;
    const Graph &graph_;
public:
    SelfComplementStat(const Graph &graph) :
            graph_(graph) {
    }

    virtual ~SelfComplementStat() {
    }

    virtual void Count() {
        size_t sc_number = 0;
        for (auto iterator = graph_.ConstEdgeBegin(); !iterator.IsEnd();
             ++iterator)
            if (graph_.conjugate(*iterator) == (*iterator))
                sc_number++;
        //    INFO("Self-complement count failed!!! ");
        INFO("Self-complement count=" << sc_number);
    }
};
}
}
