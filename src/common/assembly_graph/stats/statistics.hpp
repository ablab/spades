//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/paths/mapping_path.hpp"
#include "math/xmath.h"
#include "configs/config_struct.hpp"
#include "utils/logger/logger.hpp"
#include "utils/stl_utils.hpp"

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
    //FIXME: get rid of raw pointers
    std::vector<AbstractStatCounter *> to_count_;
public:
    StatList(const std::vector<AbstractStatCounter *> &to_count = {}) :
            to_count_(to_count) {
    }

    virtual ~StatList() {
    }

    void AddStat(AbstractStatCounter *new_stat) {
        to_count_.push_back(new_stat);
    }

    const std::vector<AbstractStatCounter *> &stats() {
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
    typedef typename Graph::EdgeId EdgeId;    
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
        for (EdgeId e : graph_.edges()) {
            edgeNumber++;
            //      if (graph_.coverage(*iterator) > 30) {
            sum_edge_length += graph_.length(e);
            //      }
        }
        return edgeNumber;
    }

    size_t edge_length() {
        size_t sum_edge_length = 0;
        for (EdgeId e : graph_.edges()) {
            if (graph_.coverage(e) > 30) {
                sum_edge_length += graph_.length(e);
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
        const auto& path_edges1 = path1_.sequence();
        const auto& path_edges2 = path2_.sequence();
        std::set<EdgeId> colored_edges;
        colored_edges.insert(path_edges1.begin(), path_edges1.end());
        colored_edges.insert(path_edges2.begin(), path_edges2.end());
        size_t sum_length = 0;
        for (EdgeId e : graph_.edges()) {
            edge_count++;
            if (colored_edges.count(e) == 0) {
                black_count++;
                sum_length += graph_.length(e);
            }
        }
        if (edge_count > 0) {
            size_t total_genome_size = 0;
            for (const auto &chr: cfg::get().ds.reference_genome)
                total_genome_size += 2*chr.size();
            INFO("Error edges count: " << black_count << " which is " <<
                 100.0 * (double) black_count / (double) edge_count << "% of all edges");
            INFO("Total length of all black edges: " << sum_length << ". While double genome length is " <<
                 total_genome_size);
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
    NStat(const Graph &graph, const Path<EdgeId> &path, size_t perc = 50) :
            graph_(graph), path_(path), perc_(perc) {
    }

    virtual ~NStat() {
    }

    virtual void Count() {
        std::vector<size_t> lengths;
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
    std::set<EdgeId> black_edges_;
    std::vector<size_t> lengths;
public:
    IsolatedEdgesStat(const Graph &graph, const Path<EdgeId> &path1, const Path<EdgeId> &path2) :
            graph_(graph) {
        for (auto e : graph_.edges())        
            black_edges_.insert(e);

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
        for (EdgeId edge : graph_.edges()) {
            if (graph_.IsDeadEnd(graph_.EdgeEnd(edge))
                && graph_.IsDeadStart(graph_.EdgeStart(edge))
                && black_edges_.count(edge) == 0) {
                lengths.push_back(graph_.length(edge));
            }
        }
        INFO("Isolated not black edges: " << lengths.size());
        WriteLengths(cfg::get().output_dir, "isolated_edges.txt");
    }

    void WriteLengths(const std::filesystem::path &folder_name, const std::filesystem::path &file_name) {
        std::ofstream os;
        os.open(folder_name / file_name);
        WriteLengths(os);
        os.close();
    }

    void WriteLengths(std::ostream &os) {
        sort(lengths.begin(), lengths.end());
        for (size_t i = 0; i < lengths.size(); i++) {
            os << lengths[i] << std::endl;
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
        for (EdgeId e : graph_.edges())
            if (graph_.conjugate(e) == (e))
                sc_number++;
        //    INFO("Self-complement count failed!!! ");
        INFO("Self-complement count=" << sc_number);
    }
};
}
}
