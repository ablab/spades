//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

namespace dipspades {

typedef vector<Sequence>::iterator cons_regions_iterator;

class ConservativeRegionStorage{
    vector<Sequence> cons_regions_;
    vector<Sequence> poss_cons_regions_;

public:
    void AddConservativeRegion(Sequence seq){
        cons_regions_.push_back(seq);
    }

    void AddPossiblyConservativeRegion(Sequence seq){
        poss_cons_regions_.push_back(seq);
    }

    cons_regions_iterator cons_regions_begin(){
        return cons_regions_.begin();
    }

    cons_regions_iterator cons_regions_end(){
        return cons_regions_.end();
    }

    cons_regions_iterator poss_cons_regions_begin(){
        return poss_cons_regions_.begin();
    }

    cons_regions_iterator poss_cons_regions_end(){
        return poss_cons_regions_.end();
    }
};

}
