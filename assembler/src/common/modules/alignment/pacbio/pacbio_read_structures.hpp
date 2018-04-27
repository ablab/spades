//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "modules/alignment/sequence_mapper.hpp"
#include "assembly_graph/core/graph.hpp"

#include <algorithm>
#include <sstream>
#include <map>
#include <set>

namespace pacbio {
typedef omnigraph::GapDescription<debruijn_graph::Graph> GapDescription;

struct MappingInstance {
    //both positions g_.k() based
    int edge_position;
    int read_position;
    //Now quality is the same with multiplicity, so best quality is 1,
    int quality;
    
    MappingInstance(int edge_position, int read_position, int quality)
            : edge_position(edge_position), read_position(read_position), quality(quality) {}

    bool IsUnique() const {
        return (quality == 1);
    }

    std::string str() {
        std::stringstream s;
        s << "E: " << edge_position << " R: " << read_position << " Q: " << quality;
        return s.str();
    }

//Less by EDGE position
    bool operator <(MappingInstance const& b) const {
        if (edge_position < b.edge_position || (edge_position == b.edge_position && read_position < b.read_position))
            return true;
        else
            return false;
    }
private:
    DECL_LOGGER("MappingInstance");
};

//Less by READ position
struct ReadPositionComparator {
    bool operator ()(MappingInstance const& a, MappingInstance const& b) const {
        return (a.read_position < b.read_position || (a.read_position == b.read_position && a.edge_position < b.edge_position));
    }
};

template<class Graph>
struct KmerCluster {
    typedef typename Graph::EdgeId EdgeId;
    int last_trustable_index;
    int first_trustable_index;
    size_t average_read_position;
    size_t average_edge_position;
    EdgeId edgeId;
    std::vector<MappingInstance> sorted_positions;
    int size;

    KmerCluster(EdgeId e, size_t edge_start_pos, size_t edge_end_pos, size_t read_start_pos, size_t read_end_pos) {
        last_trustable_index = 1;
        first_trustable_index = 0;
        sorted_positions.push_back(MappingInstance((int)edge_start_pos, (int)read_start_pos, 1));
        sorted_positions.push_back(MappingInstance((int)edge_end_pos, (int)read_end_pos, 1));
        VERIFY_MSG(edge_start_pos < edge_end_pos, "range size should be positive");
        size = int (edge_end_pos - edge_start_pos);
        average_read_position = (read_start_pos + read_end_pos)/2;
        average_edge_position = (edge_start_pos + edge_end_pos)/2;
        edgeId = e;
    }

    bool operator <(const KmerCluster & b) const {
        return (average_read_position < b.average_read_position || (average_read_position == b.average_read_position && edgeId < b.edgeId) ||
                (average_read_position == b.average_read_position && edgeId == b.edgeId && sorted_positions < b.sorted_positions));
    }

    bool CanFollow(const KmerCluster &b) const {
        return (b.sorted_positions[b.last_trustable_index].read_position < sorted_positions[first_trustable_index].read_position);
    }

    std::string str(const Graph &g) const{
        std::stringstream s;
        s << "Edge: " << g.int_id(edgeId) << " on edge: " << sorted_positions[first_trustable_index].edge_position<<
                " - "  << sorted_positions[last_trustable_index].edge_position<< ";on read: "
        << sorted_positions[first_trustable_index].read_position<< " - "
        << sorted_positions[last_trustable_index].read_position<< ";size "<< size;

        return s.str();
    }
private:
    DECL_LOGGER("KmerCluster");
};

class StatsCounter {
public:
    std::map<size_t,size_t> path_len_in_edges;
    std::vector<size_t> subreads_length;
    size_t total_len ;
    size_t reads_with_conjugate;
    size_t subreads_count;
    std::map<size_t, size_t> seeds_percentage;
    StatsCounter() {
        total_len = 0;
        reads_with_conjugate = 0;
    }

    void AddStorage(StatsCounter &other) {
        total_len += other.total_len;
        reads_with_conjugate += other.reads_with_conjugate;
        for (auto iter = other.subreads_length.begin(); iter != other.subreads_length.end(); ++iter) {
            subreads_length.push_back(*iter);
        }
        
        for (auto iter = other.path_len_in_edges.begin(); iter != other.path_len_in_edges.end(); ++iter){
            auto j_iter = iter;
            if (( j_iter = path_len_in_edges.find(iter->first)) == other.path_len_in_edges.end()){
                path_len_in_edges.insert(make_pair(iter->first, iter->second));
            } else {
                path_len_in_edges[j_iter->first] += iter->second;
            }
        }
        for (auto iter = other.seeds_percentage.begin(); iter != other.seeds_percentage.end(); ++iter){
            auto j_iter = iter;
            if (( j_iter = seeds_percentage.find(iter->first)) == other.seeds_percentage.end()){
                seeds_percentage.insert(make_pair(iter->first, iter->second));
            } else {
                seeds_percentage[j_iter->first] += iter->second;
            }
        }
    }

    void Report() const {
        size_t total = 0;
        for (auto iter = seeds_percentage.begin(); iter != seeds_percentage.end(); ++iter) {
            total += iter->second;
        }
        size_t cur = 0;
        size_t percentage = 0;
        for (auto iter = seeds_percentage.begin(); iter != seeds_percentage.end(); ++iter){
            cur += iter->second;
            percentage = iter->first;
            if (cur * 2 > total) break;
        }
        INFO("Median fraction of present seeds in maximal alignmnent among reads aligned to the graph: " << double(percentage) * 0.001);
    }

private:
    DECL_LOGGER("StatsCounter");
};

}
