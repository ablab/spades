//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "redundant_contig_remover.hpp"
#include "equal_path_deletion_correction.hpp"

using namespace debruijn_graph;

namespace dipspades {

class IterativeLoopCorrector : public AbstractContigCorrector{

    size_t k_value_;

    VertexPathIndex &index_;
    size_t max_loop_len_;
    size_t min_lcs_length_;
    size_t max_tail_length_;
    CorrectionResult res;

public:
    IterativeLoopCorrector(Graph &g, size_t k_value, VertexPathIndex &index, size_t max_loop_len,
            size_t min_lcs_length, size_t max_tail_length) :
        AbstractContigCorrector(g), k_value_(k_value), index_(index),
        max_loop_len_(max_loop_len), min_lcs_length_(min_lcs_length),
        max_tail_length_(max_tail_length) {
    }

    ContigStoragePtr Correct(ContigStoragePtr contigs) {
        {
            INFO("Equal path remover starts");
            index_.Initialize(contigs);
            EqualPathDeletionCorrector equal_path_remover(g_, index_);
            contigs = equal_path_remover.Correct(contigs);
            res.redundancy_map = equal_path_remover.Result().redundancy_map;
            index_.Clear();
            INFO(ToString(contigs->Size()) + " contigs will be used further");
        }

        INFO("Iterative loop corrector starts");
        {
            INFO("Only exact match iteration with parameters:");
            INFO("\tMaximal loop length - " + ToString(max_loop_len_));
            INFO("\tMinimal lcs length - " + ToString(min_lcs_length_));
            INFO("\tMaximal tail length - 0");

            index_.Initialize(contigs);
            LoopBulgeDeletionCorrector loop_corr(g_, k_value_,
                    max_loop_len_, 0, min_lcs_length_, index_);
            contigs = loop_corr.Correct(contigs);
            auto old_map = res.redundancy_map;
            auto new_map = loop_corr.Results().redundancy_map;
            RedundancyMapMerger<size_t> map_merger;
            res.redundancy_map = map_merger.MergeTwoMaps(old_map, new_map);
            index_.Clear();
            INFO(ToString(contigs->Size()) + " contigs will be used further");
        }

        {
            INFO("Tails allowing match iteration with parameters:");
            INFO("\tMaximal loop length - " + ToString(max_loop_len_));
            INFO("\tMinimal lcs length - " + ToString(min_lcs_length_));
            INFO("\tMaximal tail length - " + ToString(max_tail_length_));
            index_.Initialize(contigs);
            LoopBulgeDeletionCorrector loop_corr(g_, k_value_,
                    max_loop_len_, max_tail_length_, min_lcs_length_, index_);
            contigs = loop_corr.Correct(contigs);
            auto old_map = res.redundancy_map;
            auto new_map = loop_corr.Results().redundancy_map;
            RedundancyMapMerger<size_t> map_merger;
            res.redundancy_map = map_merger.MergeTwoMaps(old_map, new_map);
            index_.Clear();
            INFO(ToString(contigs->Size()) + " contigs will be used further");
        }
        INFO("Iterative loop corrector ends");
        return contigs;
    }

    MappingContigPtr Correct(MappingContigPtr contig){
        return contig;
    }

    CorrectionResult Results(){
        return res;
    }
};

}
