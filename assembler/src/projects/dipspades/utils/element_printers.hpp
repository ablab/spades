//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <set>

using namespace std;
using namespace debruijn_graph;

namespace dipspades {

void PrintSimplePath(ostream &out, Graph &g, vector<EdgeId> path){
    for(auto e = path.begin(); e != path.end(); e++)
        out << g.int_id(*e) << " ";
    out << endl;
}

void PrintSimplePathWithVertices(ostream &out, Graph &g, vector<EdgeId> &path){
    for(auto e = path.begin(); e != path.end(); e++)
        out << g.int_id(*e) << " (" << g.length(*e) << "), " << g.int_id(g.EdgeStart(*e)) << " - " <<
        g.int_id(g.EdgeEnd(*e)) << ". ";
    out << endl;
}

string SimplePathWithVerticesToString(const Graph &g, vector<EdgeId> path){
    stringstream ss;
    for(auto e = path.begin(); e != path.end(); e++)
        ss << g.int_id(*e) << " (" << g.length(*e) << "), " << g.int_id(g.EdgeStart(*e)) << " - " <<
        g.int_id(g.EdgeEnd(*e)) << ". ";
    return ss.str();
}

string MappingPathToString(Graph &g, MappingPath<EdgeId> path){
    stringstream ss;
    for(size_t i = 0; i < path.size(); i++){
        Range init = path[i].second.initial_range, mapp = path[i].second.mapped_range;
        ss << "Edge - " << g.str(path[i].first) <<  " (" << g.length(path[i].first) << ") . Init range - " << init.start_pos <<
                " - " << init.end_pos << ". Mapp range - " << mapp.start_pos << " - " <<
                mapp.end_pos << ". ";
    }
    return ss.str();
}

template<class T>
void PrintSet(ostream &out, set<T> set_elem){
    for(auto e = set_elem.begin(); e != set_elem.end(); e++)
        out << *e << " ";
    out << endl;
}

template<class T>
string SetToString(set<T> set_elem){
    stringstream ss;
    for(auto e = set_elem.begin(); e != set_elem.end(); e++)
        ss << *e << " ";
    return ss.str();
}

template<class T>
void PrintVector(ostream &out, vector<T> vect_elem){
    for(auto e = vect_elem.begin(); e != vect_elem.end(); e++)
        out << *e << " ";
    out << endl;
}

template<class T>
string VectorToString(vector<T> vect_elem){
    stringstream ss;
    for(auto e = vect_elem.begin(); e != vect_elem.end(); e++)
        ss << *e << " ";
    return ss.str();
}

string VerticesVectorToString(Graph &g, vector<VertexId> vertices){
    stringstream ss;
    for(auto it = vertices.begin(); it != vertices.end(); it++)
        ss << g.str(*it) << " ";
    return ss.str();
}

void PrintEdgeWithVertices(ostream &out, Graph &g, EdgeId edge){
    out << "Edge - " << g.int_id(edge) << ". Start vertex - " << g.int_id(g.EdgeStart(edge)) <<
            ". End vertex - " << g.int_id(g.EdgeEnd(edge)) << endl;
}

void PrintEdgeWithLength(ostream &out, Graph &g, EdgeId edge){
    out << "Edge - " << g.int_id(edge) << " with length - " << g.length(edge) << endl;
}

void PrintVectorOfVertices(ostream &out, Graph &g, vector<VertexId> vect){
    for(auto v = vect.begin(); v != vect.end(); v++)
        out << g.int_id(*v) << " ";
    out << endl;
}

string MappingRangeToString(MappingRange mr){
    stringstream ss;
    ss << "Init: " << mr.initial_range.start_pos << " " << mr.initial_range.end_pos
            << ". Map: " << mr.mapped_range.start_pos << " " << mr.mapped_range.end_pos << endl;
    return ss.str();
}

}
