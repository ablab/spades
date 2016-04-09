//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "abstract_contig_corrector.hpp"

using namespace debruijn_graph;

namespace dipspades {

class EqualPathDeletionCorrector : public AbstractContigCorrector{
    VertexPathIndex &path_index_;
    CorrectionResult res_;

    void InitializeMap(ContigStoragePtr contigs){
        for(size_t i = 0; i < contigs->Size(); i++){
            size_t id = (*contigs)[i]->id();
            res_.redundancy_map.AddNewKey(id);
        }
    }

public:

    EqualPathDeletionCorrector(Graph &g, VertexPathIndex &path_index) : AbstractContigCorrector(g),
    path_index_(path_index){ }

    ContigStoragePtr Correct(ContigStoragePtr contigs)    {

        INFO("Computing redundant equal contigs starts");

        InitializeMap(contigs);
        set<size_t> ids_for_deletion;
        for(size_t i = 0; i < contigs->Size() - 1; i++){
            size_t id1 = (*contigs)[i]->id();
            size_t rc_id1 = (*contigs)[i]->rc_id();
            if(ids_for_deletion.find(id1) == ids_for_deletion.end() &&
                    ids_for_deletion.find(rc_id1) == ids_for_deletion.end()){
                auto path1 = (*contigs)[i]->path_seq();
                auto contigs_for_processing = path_index_.GetPathsIntersectedWith(path1);
                for(auto it = contigs_for_processing.begin(); it != contigs_for_processing.end(); it++){
                    size_t j = *it;
                    size_t id2 = (*contigs)[j]->id();
                    size_t rc_id2 = (*contigs)[j]->rc_id();
                    if(ids_for_deletion.find(id2) == ids_for_deletion.end() &&
                            ids_for_deletion.find(rc_id2) == ids_for_deletion.end() && j > i){
                        auto path2 = (*contigs)[j]->path_seq();
                        if(ArePathEqual(path1, path2)){
                            size_t id2 = (*contigs)[j]->id();
                            ids_for_deletion.insert(id2);
                            ids_for_deletion.insert(rc_id2);
                            res_.redundancy_map.AddNewPair(id1, id2);
                            res_.redundancy_map.AddNewPair(rc_id1, rc_id2);
                        }
                    }
                }
            }
        }
        RedundancyMapCondenser<size_t> condenser;
        res_.redundancy_map = condenser.Condense(res_.redundancy_map);
        INFO(ToString(ids_for_deletion.size()) + " contigs from " << contigs->Size() << " are redundant");
        contigs->DeleteByIDs(ids_for_deletion);

        INFO("Computing redundant equal contigs ends");

        return contigs;
    }

    MappingContigPtr Correct(MappingContigPtr contig){
        return contig;
    }

    CorrectionResult Result(){
        return res_;
    }
};

}
