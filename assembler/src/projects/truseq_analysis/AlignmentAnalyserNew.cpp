//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

//
// Created by anton on 5/15/15.
//

#include "utils/standard_base.hpp"
#include "assembly_graph/dijkstra/dijkstra_helper.hpp"
#include "AlignmentAnalyserNew.hpp"

namespace alignment_analysis {
    using omnigraph::Range;

    size_t AlignmentAnalyserNew::StepBack(const vector<ConsistentMapping> &path) const {
        size_t cur_step = 0;
        size_t cur = path.size() - 1;
        while(cur > 0 && cur_step + path[cur].size() < step_ && path[cur - 1].CheckConnect(path[cur])) {
            cur_step += path[cur].size();
            cur--;
        }
        return cur;
    }

    vector <ConsistentMapping> AlignmentAnalyserNew::Analyse(const omnigraph::MappingPath<EdgeId> &path) const {
        TRACE("");
        TRACE("Analysis of path of length " << path.size() << ": " << path);
        vector <ConsistentMapping> result;
        TRACE("Adding " << path[0]);
        result.push_back(ConsistentMapping(graph_, path[0].first, path[0].second));
        for(size_t i = 1; i < path.size(); i++) {
            TRACE("Attempting to add new part");
            ConsistentMapping &mapping = result.back();
            if(mapping.CheckConnect(path[i].first, path[i].second)) {
                TRACE("Adding #" << result.size() << ": " << path[i]);
                result.push_back(ConsistentMapping(graph_, path[i].first, path[i].second));
            } else if (mapping.EndEdge() == path[i].first && mapping.Back().second.end_pos <= path[i].second.mapped_range.start_pos) {
                Range initial(mapping.GetInitialRange().end_pos, path[i].second.initial_range.start_pos);
                Range mapped(mapping.Back().second.end_pos, path[i].second.mapped_range.start_pos);
                result.push_back(ConsistentMapping(graph_, path[i].first, omnigraph::MappingRange(initial, mapped)));
                result.push_back(ConsistentMapping(graph_, path[i].first, path[i].second));
            } else {
                TRACE("Could not add " << path[i]);
                size_t pos = StepBack(result);
                VertexId start = result[pos].EndVertex();
                TRACE("Start vertex: " << start);
                omnigraph::DijkstraHelper<Graph>::BoundedDijkstra d = omnigraph::DijkstraHelper<Graph>::CreateBoundedDijkstra(graph_, 3000 + graph_.k(), 1000);
                d.Run(start);
                size_t best = i;
                for (size_t j = i, cur_step = 0; j < path.size() && cur_step < step_; j++) {
                    TRACE("Checking candidate #" << j << ": " << path[j]);
                    if (d.DistanceCounted(graph_.EdgeStart(path[j].first))){
                        best = j;
                        break;
                    }
                    cur_step += path[j].second.mapped_range.size();
                }
                if(best < path.size() && d.DistanceCounted(graph_.EdgeStart(path[best].first))) {
                    //todo: make better cutting
                    this->Cut(result, start);
                    vector<EdgeId> detour_path = d.GetShortestPathTo(graph_.EdgeStart(path[best].first));
                    vector<EdgeRange> detour_mapping;
                    EdgeRange er = result.back().GetMappedPath().back();
                    if (er.second.end_pos != graph_.length(er.first)) {
                        detour_mapping.push_back(EdgeRange(er.first, Range(er.second.end_pos, graph_.length(er.first))));
                    }
                    for(auto it = detour_path.begin(); it != detour_path.end(); ++it) {
                        detour_mapping.push_back(EdgeRange(*it, Range(0, graph_.length(*it))));
                    }
                    if (path[best].second.mapped_range.start_pos != 0) {
                        detour_mapping.push_back(EdgeRange(path[best].first, Range(0, path[best].second.mapped_range.start_pos)));
                    }
                    if(detour_mapping.size() > 0) {
                        Range r(result.back().GetInitialRange().end_pos, path[best].second.initial_range.start_pos);
                        ConsistentMapping detour(graph_, r, detour_mapping);
                        TRACE("Adding #" << result.size() << ": " << detour);
                        result.push_back(detour);
                    } else {
                        TRACE("Empty detour path");
                    }
                    TRACE("Adding #" << result.size() << ": " << path[best]);
                    result.push_back(ConsistentMapping(graph_, path[best].first, path[best].second));
                    i = best;
                } else {
                    TRACE("Adding #" << result.size() << ": " << path[i]);
                    result.push_back(ConsistentMapping(graph_, path[i].first, path[i].second));
                }
            }
        }
        return result;
    }

    void AlignmentAnalyserNew::Cut(vector<ConsistentMapping> &path, VertexId start) const {
        while(path.back().EndVertex() != start) {
            path.pop_back();
        }
    }
}
