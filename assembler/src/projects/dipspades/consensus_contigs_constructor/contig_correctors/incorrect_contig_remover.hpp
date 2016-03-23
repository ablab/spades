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

class RemoveUnconnectContigsCorrector : public AbstractContigCorrector{

public:
    RemoveUnconnectContigsCorrector(Graph &g) : AbstractContigCorrector(g){ }

    ContigStoragePtr Correct(ContigStoragePtr storage) {
        set<size_t> contigs_for_deletion;
        for(size_t i = 0; i < storage->Size(); i++){
            auto contig_path = (*storage)[i]->path_seq();
            TRACE((*storage)[i]->id() << " contig");
            TRACE("Path: " << SimplePathWithVerticesToString(g_, contig_path));
            if(!IsPathConnected(g_, contig_path)){
                contigs_for_deletion.insert((*storage)[i]->id());
            }
        }
        INFO(ToString(contigs_for_deletion.size()) +  " contigs from " <<
                storage->Size() << " were deleted");
        storage->DeleteByIDs(contigs_for_deletion);
        return storage;
    }

    MappingContigPtr Correct(MappingContigPtr contig){
        return contig;
    }

};

}
