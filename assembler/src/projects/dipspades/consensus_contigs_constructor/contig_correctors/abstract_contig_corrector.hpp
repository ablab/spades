//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "../mapping_contigs_storage.hpp"
#include "../../utils/lcs_utils.hpp"
#include "../../utils/path_routines.hpp"
#include "../../utils/path_index.hpp"
#include "../../utils/bulge_utils.hpp"
#include "../../utils/redundancy_map.hpp"

using namespace debruijn_graph;

namespace dipspades {

struct CorrectionResult{
    OverlapGraph g;
    RedundancyMap<size_t> redundancy_map;
};

//--------------------------------------------------------------------------
class AbstractContigCorrector{
protected:
    Graph& g_;
public:
    AbstractContigCorrector(Graph& g) : g_(g) {

    }
    virtual ContigStoragePtr Correct(ContigStoragePtr storage) { return storage; }
    virtual MappingContigPtr Correct(MappingContigPtr contig) { return contig; }
    virtual ~AbstractContigCorrector(){}
    virtual CorrectionResult Results(){
        CorrectionResult res;
        return res;
    }
};

}
