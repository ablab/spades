//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once
#include "assembly_graph/core/graph.hpp"
#include "pipeline/graph_pack.hpp"

using namespace debruijn_graph;

namespace dipspades {

bool AreEdgesConnected(Graph &g, EdgeId e1, EdgeId e2){
    return g.EdgeEnd(e1) == g.EdgeStart(e2);
}

bool IsPathConnected(Graph &g, vector<EdgeId> path){
    if(path.size() <= 1)
        return true;
    for(size_t i = 0; i < path.size() - 1; i++){
        EdgeId e1 = path[i];
        EdgeId e2 = path[i + 1];
        if(!AreEdgesConnected(g, e1, e2)){
            return false;
        }
    }
    return true;
}

bool PathContainsLoop(Graph &g, vector<EdgeId> p){
    if(p.size() == 0)
        return false;

    set<VertexId> pathv;
    pathv.insert(g.EdgeStart(p[0]));

    for(auto e = p.begin(); e != p.end(); e++){
        VertexId end = g.EdgeEnd(*e);
        if(pathv.find(end) == pathv.end())
            pathv.insert(end);
        else
            return true;
    }

    return false;
}

vector<VertexId> get_list_of_vertices_in_path(Graph &g, vector<EdgeId> path){
    vector<VertexId> list;
    if(path.size() == 0)
        return list;

    for(size_t i = 0; i < path.size(); i++)
        list.push_back(g.EdgeStart(path[i]));

    list.push_back(g.EdgeEnd(path[path.size() - 1]));

    return list;
}

bool is_1st_edge_not_later_2nd_edge(vector<EdgeId> path, EdgeId first_edge,
        EdgeId second_edge){
    bool first_found = false;
    for(auto e = path.begin(); e != path.end(); e++){
        if(*e == first_edge)
            first_found = true;
        if(*e == second_edge && !first_found)
            return false;
        if(*e == second_edge && first_found)
            return true;
    }

    return false;
}

bool is_1st_edge_not_later_2nd_edge(vector<EdgeId> path, EdgeId first_edge,
        EdgeId second_edge, int ind1, int ind2){
    bool first_found = false;
    //for(auto e = path.begin(); e != path.end(); e++){
    for(int i = ind1; i <= ind2; i++){
        EdgeId e = path[i];
        if(e == first_edge)
            first_found = true;
        if(e == second_edge && !first_found)
            return false;
        if(e == second_edge && first_found)
            return true;
    }

    return false;
}

int get_index_of_edge(vector<EdgeId> path, EdgeId edge){
    for(size_t i = 0; i < path.size(); i++)
        if(path[i] == edge)
            return int(i);
    return -1;
}

EdgeId GetEdgeById(conj_graph_pack & gp, size_t id){
    for(auto e = gp.g.SmartEdgeBegin(); !e.IsEnd(); ++e)
        if(gp.g.int_id(*e) == id)
            return *e;
    return EdgeId(0);
}

bool IsPathRegionCorrect(pair<size_t, size_t> region, size_t path_size){
    return region.first < path_size && region.second < path_size;
}

size_t GetLengthOfPathRegion(Graph &g, vector<EdgeId> path, pair<size_t, size_t> region){
    VERIFY(IsPathRegionCorrect(region, path.size()));
    size_t region_length = 0;
    for(size_t i = region.first; i <= region.second; i++)
        region_length += g.length(path[i]);
    return region_length;
}

size_t GetPathLength(Graph &g, vector<EdgeId> path){
    if(path.size() == 0)
        return 0;
    return GetLengthOfPathRegion(g, path, pair<size_t, size_t>(0, path.size() - 1));
}

Sequence GetSequenceOfPathRegion(Graph &g, size_t k_value, vector<EdgeId> path,
        pair<size_t, size_t> region){
    VERIFY(IsPathRegionCorrect(region, path.size()));

    if(region.first > region.second)
        return Sequence();

    EdgeId cur_edge = path[region.first];
    Sequence seq = g.EdgeNucls(cur_edge);

    for(auto i = region.first + 1; i <= region.second; ++i){
        Sequence next_seq = g.EdgeNucls(path[i]);
        seq = seq + next_seq.Subseq(k_value, next_seq.size());
    }

    return seq;
}

Sequence GetSequenceByPath(Graph &g, size_t k_value, const vector<EdgeId> path){
    if(path.size() == 0)
        return Sequence();
    return GetSequenceOfPathRegion(g, k_value, path, pair<size_t, size_t>(0, path.size() - 1));
}

Sequence GetSequenceByPath(conj_graph_pack &gp, const vector<EdgeId> path){
    if(path.size() == 0)
        return Sequence();
    return GetSequenceOfPathRegion(gp.g, gp.k_value, path, pair<int, int>(0, path.size() - 1));
}

vector<EdgeId> GetRCToPathSeq(Graph &g, vector<EdgeId> path){
    vector<EdgeId> rc_path;
    for(auto e = path.begin(); e != path.end(); e++){
        rc_path.insert(rc_path.begin(), g.conjugate(*e));
    }
    return rc_path;
}

MappingPath<EdgeId> GetRCToMappingPath(Graph &g, MappingPath<EdgeId> map_path, size_t seq_size){
    vector<EdgeId> rc_path_seq;
    vector<MappingRange> rc_map_ranges;
    for(size_t i = 0; i < map_path.size(); i++){
        // computing edges sequence
        EdgeId cur_edge = map_path[i].first;
        rc_path_seq.insert(rc_path_seq.begin(), g.conjugate(cur_edge));

        // computing initial ranges
        Range init_range = map_path[i].second.initial_range;
        Range rc_init_range(seq_size - init_range.end_pos, seq_size - init_range.start_pos);

        // computing mapped ranges
        size_t edge_length = g.length(cur_edge);
        Range map_range = map_path[i].second.mapped_range;
        Range rc_map_range(edge_length - map_range.end_pos, edge_length - map_range.start_pos);

        rc_map_ranges.insert(rc_map_ranges.begin(), MappingRange(rc_init_range, rc_map_range));
    }

    return MappingPath<EdgeId>(rc_path_seq, rc_map_ranges);
}

bool ArePathEqual(vector<EdgeId> path1, vector<EdgeId> path2){
    if(path1.size() != path2.size())
        return false;

    for(size_t i = 0; i < path1.size(); i++)
        if(path1[i] != path2[i])
            return false;

    return true;
}

bool PathsShareEdge(vector<EdgeId> path1, vector<EdgeId> path2){
    for(auto it1 = path1.begin(); it1 != path1.end(); it1++)
        for(auto it2 = path2.begin(); it2 != path2.end(); it2++)
            if(*it1 == *it2)
                return true;
    return false;
}

vector<EdgeId> CutSubpathByRegion(vector<EdgeId> path, pair<size_t, size_t> region){
    VERIFY(IsPathRegionCorrect(region, path.size()));
    vector<EdgeId> subpath;
    for(size_t i = region.first; i <= region.second; i++)
        subpath.push_back(path[i]);
    return subpath;
}

bool IsEdgeRelated(Graph &g, EdgeId edge){
    return g.RelatedVertices(g.EdgeStart(edge), g.EdgeEnd(edge));
}

bool IsEdgeLoop(Graph &g, EdgeId edge){
    return g.EdgeStart(edge) == g.EdgeEnd(edge);
}

bool VertexAdjacentRelatedEdges(Graph &g, VertexId vertex){
    auto in_edges = g.IncomingEdges(vertex);
    for(auto it = in_edges.begin(); it != in_edges.end(); it++)
        if(IsEdgeRelated(g, *it))
            return true;
    auto out_edges = g.OutgoingEdges(vertex);
    for(auto it = out_edges.begin(); it != out_edges.end(); it++)
        if(IsEdgeRelated(g, *it))
            return true;
    return false;
}

bool PathAdjacentRelatedEdges(Graph &g, vector<EdgeId> path, bool check_start = false,
        bool check_end = false){
    for(auto e = path.begin(); e != path.end() - 1; e++)
        if(VertexAdjacentRelatedEdges(g, g.EdgeEnd(*e)))
            return true;
    if(path.size() != 0) {
        if(check_start)
            if(VertexAdjacentRelatedEdges(g, g.EdgeStart(path[0])))
                return true;
        if(check_end)
            if(VertexAdjacentRelatedEdges(g, g.EdgeEnd(path[path.size() - 1])))
                return true;
    }
    return false;
}

vector<size_t> CalculatePathPartLens(Graph &g, vector<EdgeId> path){
    vector<size_t> lens;
    size_t prev_len = 0;
    for(auto e = path.begin(); e != path.end(); e++){
        lens.push_back(prev_len + g.length(*e));
        prev_len += g.length(*e);
    }
    return lens;
}

/*
void detect_loop_length(ostream &out, Graph& g, ContigStorage* stor){

    for(size_t i = 0; i < stor->Size(); i++){
        vector<EdgeId> path = (*stor)[i]->PathSeq();

        vector<VertexId> vert =  get_list_of_vertices_in_path(g, path);

        for(size_t j = 0; j < vert.size() - 1; j++){
            for(size_t k = j + 1; k < vert.size(); k++){
                if(vert[j] == vert[k]){
                    size_t ind1 = j, ind2 = (k == path.size()) ? (k - 1) : k;
                    size_t loop_length = 0;

                    for(size_t l = ind1; l <= ind2; l++)
                        loop_length += g.length(path[l]);
                    out << loop_length << endl;
                }
            }
        }
    }
}
*/

}
