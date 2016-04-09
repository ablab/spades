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

class SameEdgeDeletionCorrector : public AbstractContigCorrector{

    MappingRange WideMappingRange(MappingRange old_range, MappingRange new_range){

        size_t initial_range_start_pos = min<size_t>(old_range.initial_range.start_pos, new_range.initial_range.start_pos);
        size_t initial_range_end_pos = max<size_t>(old_range.initial_range.end_pos, new_range.initial_range.end_pos);
        size_t mapped_range_start_pos = min<size_t>(old_range.mapped_range.start_pos, new_range.mapped_range.start_pos);
        size_t mapped_range_end_pos = max<size_t>(old_range.mapped_range.end_pos, new_range.mapped_range.end_pos);

        Range init(initial_range_start_pos, initial_range_end_pos), mapp(mapped_range_start_pos, mapped_range_end_pos);
        MappingRange res(init, mapp);
        return res;
    }

public:
    SameEdgeDeletionCorrector(Graph &g) : AbstractContigCorrector(g) {
    }

    ContigStoragePtr Correct(ContigStoragePtr contigs) {
        for(size_t i = 0; i < contigs->Size(); i++)
            contigs->ReplaceContig(Correct((*contigs)[i]), i);
        TRACE(contigs->Size() << " contigs from " << contigs->Size() << " were corrected");
        return contigs;
    }

    MappingContigPtr Correct(MappingContigPtr contig){
        MappingPath<EdgeId> map_path = contig->mapping_path();

        if(map_path.size() <= 0)
            return contig;

        vector<EdgeId> new_path;
        vector<MappingRange> new_ranges;
        EdgeId cur_edge = map_path[0].first;
        new_path.push_back(cur_edge);
        new_ranges.push_back(map_path[0].second);

        for (size_t i = 1; i < map_path.size(); i++) {
            EdgeId e = map_path[i].first;
            if (e != cur_edge) {
                cur_edge = e;
                new_path.push_back(e);
                new_ranges.push_back(map_path[i].second);
            }
            else {
                new_ranges[new_ranges.size() - 1] = WideMappingRange(
                        new_ranges[new_ranges.size() - 1], map_path[i].second);
            }
        }

        MappingPath<EdgeId> new_map_path(new_path, new_ranges);
        return MappingContigPtr(new ReplacedPathMappingContig(contig, new_map_path));
    }
};

}
