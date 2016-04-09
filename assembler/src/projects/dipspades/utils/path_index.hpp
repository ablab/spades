//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "../consensus_contigs_constructor/mapping_contigs_storage.hpp"

namespace dipspades {

class VertexPathIndex{
    Graph &g_;

    map<VertexId, set<size_t> > index_;

    void AddNewPair(VertexId v, size_t path_index){
        index_[v].insert(path_index);
    }

    set<size_t> JoinTwoSets(set<size_t> set1, set<size_t> set2){
        for(auto it = set2.begin(); it != set2.end(); it++)
            set1.insert(*it);
        return set1;
    }

public:
    VertexPathIndex(Graph &g) : g_(g) {
    }

    void Initialize(ContigStoragePtr storage){
        INFO("Initialization of vertex-paths index starts");
        for(size_t i = 0; i < storage->Size(); i++){
            auto path = (*storage)[i]->path_seq();
            if(path.size() > 0){
                VertexId start_vertex = g_.EdgeStart(path[0]);
                AddNewPair(start_vertex, i);

                for(auto e = path.begin(); e != path.end(); e++){
                    VertexId v = g_.EdgeEnd(*e);
                    AddNewPair(v, i);
                }
            }
        }
        INFO("Initialization of vertex-paths index ends");
    }

    void Clear(){
        index_.clear();
    }

    set<size_t> GetPathsIntersectedWith(vector<EdgeId> path){
        set<size_t> res;
        if(path.size() == 0)
            return res;
        VertexId start_vertex = g_.EdgeStart(path[0]);
        res = index_[start_vertex];
        for(auto e = path.begin(); e != path.end(); e++){
            VertexId v = g_.EdgeEnd(*e);
            res = JoinTwoSets(res, index_[v]);
        }
        return res;
    }
};

}
