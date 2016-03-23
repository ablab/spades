//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

using namespace debruijn_graph;
using namespace std;

#include "bulge_utils.hpp"

namespace dipspades {

double RelativeAlignmentOfSequencesByMinimal(Sequence seq1, Sequence seq2, bool from_start = true){
    Sequence trim_seq1, trim_seq2;
    if(min<size_t>(seq1.size(), seq2.size()) == seq1.size()){
        trim_seq1 = seq1;
        if(from_start)
            trim_seq2 = seq2.Subseq(seq2.size());
        else
            trim_seq2 = seq2.Subseq(seq2.size() - seq1.size(), seq1.size());
    }
    else{
        if(from_start)
            trim_seq1 = seq1.Subseq(seq2.size());
        else
            trim_seq1 = seq1.Subseq(seq1.size() - seq2.size(), seq1.size());
        trim_seq2 = seq2;
    }
    return RelAlignmentOfSequences(trim_seq1, trim_seq2);
}

}
