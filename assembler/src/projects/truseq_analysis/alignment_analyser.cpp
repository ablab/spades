//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "utils/standard_base.hpp"
#include "alignment_analyser.hpp"

namespace alignment_analysis {
    AlignmentAnalyser::AlignmentAnalyser(const vector <io::SingleRead> &scaffolds,
                                         const vector <io::SingleRead> &genome,
                                         const Graph &graph, const Mapper &mapper) : graph_(graph),
                                                                                     mapper_(mapper),
                                                                                     scaffolds_(scaffolds),
                                                                                     genome_(genome) {
    }

    string AlignmentAnalyser::str(const EdgeRange &er) const {
        stringstream result;
        result << "[" << er.first.int_id() << " len: " << er.first->length(55) << " from: " <<
        graph_.EdgeStart(er.first) << " to: " << graph_.EdgeEnd(er.first) << ", " << er.second << "]";
        return result.str();
    }

    using std::cout;
    using std::endl;

    vector <ConsistentMapping> AlignmentAnalyser::DetectAndMaskShortMutations(
            const vector <ConsistentMapping> &alignments) {

        if (alignments.empty()) {
            return vector<ConsistentMapping>();
        }
        DijkstraHelper<Graph>::BoundedDijkstra d = DijkstraHelper<Graph>::CreateBoundedDijkstra(graph_,
                                                                                                3000 + graph_.k(),
                                                                                                1000);
        vector <ConsistentMapping> result = {alignments.front()};
        for (size_t i = 0; i + 1 < alignments.size(); i++) {
            ConsistentMapping &prev = result.back();
            const ConsistentMapping &next = alignments[i + 1];
            const EdgeRange &back = prev.Back();
            const EdgeRange &front = next.Front();
            if (back.first == front.first) {
                vector <EdgeRange> v;
                if (back.second.end_pos < front.second.start_pos)
                    v.push_back(EdgeRange(back.first,
                                                             Range(back.second.end_pos, front.second.start_pos)));
                if (back.second.end_pos <= front.second.start_pos) {
                    prev.Join(next, v);
                    continue;
                } else {
                    std::cout << "incompatible alignments" << std::endl;
                }
            }
            VertexId gap_start = graph_.EdgeEnd(back.first);
            VertexId gap_end = graph_.EdgeStart(front.first);
            d.Run(gap_start);
            if (d.DistanceCounted(gap_end)) {
                vector <EdgeId> path = d.GetShortestPathTo(gap_end);
                int s = (int(graph_.length(back.first)) - int(back.second.end_pos)) + int(front.second.start_pos);
                for (auto it = path.begin(); it != path.end(); ++it)
                    s += graph_.length(*it);
                int diff = int(next.GetInitialRange().start_pos) - int(prev.GetInitialRange().end_pos) - s;
                this->log_ << "Found short mutation: segment [" << prev.GetInitialRange().end_pos << ", " <<
                next.GetInitialRange().start_pos <<
                "] was replaced with path of length " << s << "(difference : " << diff << ")" << endl;
                result.back().ForceJoin(next, path);
            } else {
                result.push_back(next);
                if (graph_.OutgoingEdgeCount(graph_.EdgeEnd(back.first)) == 0 &&
                    graph_.IncomingEdgeCount(graph_.EdgeStart(front.first)) == 0
                    && back.second.end_pos + 100 > graph_.length(back.first) && front.second.start_pos < 100) {
                    this->log_ << "Coverage break detected " << prev.GetInitialRange().end_pos << " till " <<
                    next.GetInitialRange().start_pos << endl;
                } else {
                    this->log_ << "Found break: from " << prev.GetInitialRange().end_pos << " till " <<
                    next.GetInitialRange().start_pos << endl;
                }
            }
        }
        return result;
    }

    vector <ConsistentMapping> AlignmentAnalyser::ExtractConsistentMappings(const MappingPath<EdgeId> &path) {
        vector <ConsistentMapping> result;
        for (size_t i = 0; i < path.size(); i++) {
            pair<EdgeId, MappingRange> m = path[i];
            ConsistentMapping mapping = ConsistentMapping(this->graph_, m.first, m.second);
            if (result.empty()) {
                result.push_back(mapping);
            } else {
                ConsistentMapping &back = result.back();
                if (back.CheckConnect(m.first, m.second.mapped_range) &&
                    back.GetInitialRange().end_pos == m.second.initial_range.start_pos) {
                    back.Join(mapping);
                } else {
                    result.push_back(mapping);
                }
            }
        }
        return result;
    }

    string AlignmentAnalyser::Analyse(const io::SingleRead &genome_part) {
        log_.str("");
        MappingPath<EdgeId> path = mapper_.MapRead(genome_part);
        stringstream result;
        log_ << "Analysis of part " << genome_part.name() << endl;
        cout << "Analysis of part " << genome_part.name() << endl;
        vector <ConsistentMapping> mapping = ExtractConsistentMappings(path);
        mapping = DetectAndMaskShortMutations(mapping);
        return log_.str();
    }
}
