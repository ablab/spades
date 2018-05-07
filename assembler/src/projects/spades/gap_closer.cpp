//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "gap_closer.hpp"
#include "assembly_graph/stats/picture_dump.hpp"
#include "modules/simplification/compressor.hpp"
#include "io/dataset_support/read_converter.hpp"
#include <stack>

namespace debruijn_graph {

class GapCloserPairedIndexFiller {
private:
    const Graph &graph_;
    const SequenceMapper<Graph> &mapper_;

    size_t CorrectLength(Path<EdgeId> path, size_t idx) const {
        size_t answer = graph_.length(path[idx]);
        if (idx == 0)
            answer -= path.start_pos();
        if (idx == path.size() - 1)
            answer -= graph_.length(path[idx]) - path.end_pos();
        return answer;
    }

    template<typename PairedRead>
    void ProcessPairedRead(omnigraph::de::PairedInfoBuffer<Graph> &paired_index,
                           const PairedRead &p_r,
                           const std::unordered_map<EdgeId, pair<EdgeId, int> > &OutTipMap,
                           const std::unordered_map<EdgeId, pair<EdgeId, int> > &InTipMap) const {
        Sequence read1 = p_r.first().sequence();
        Sequence read2 = p_r.second().sequence();

        Path<EdgeId> path1 = mapper_.MapSequence(read1).path();
        Path<EdgeId> path2 = mapper_.MapSequence(read2).path();
        for (size_t i = 0; i < path1.size(); ++i) {
            auto OutTipIter = OutTipMap.find(path1[i]);
            if (OutTipIter != OutTipMap.end()) {
                for (size_t j = 0; j < path2.size(); ++j) {
                    auto InTipIter = InTipMap.find(path2[j]);
                    if (InTipIter != InTipMap.end()) {
                        auto e1 = OutTipIter->second.first;
                        auto e2 = InTipIter->second.first;
                        //FIXME: Normalize fake points
                        auto sp = std::make_pair(e1, e2);
                        auto cp = paired_index.ConjugatePair(e1, e2);
                        auto ip = std::min(sp, cp);
                        paired_index.Add(ip.first, ip.second, omnigraph::de::RawPoint(1000000., 1.));
                    }
                }
            }
        }
    }

    void PrepareShiftMaps(std::unordered_map<EdgeId, pair<EdgeId, int> > &OutTipMap,
                          std::unordered_map<EdgeId, pair<EdgeId, int> > &InTipMap) {
        std::stack<pair<EdgeId, int>> edge_stack;
        for (auto iterator = graph_.ConstEdgeBegin(); !iterator.IsEnd();) {
            EdgeId edge = *iterator;
            if (graph_.IncomingEdgeCount(graph_.EdgeStart(edge)) == 0) {
                InTipMap.insert(std::make_pair(edge, std::make_pair(edge, 0)));
                edge_stack.push(std::make_pair(edge, 0));
                while (edge_stack.size() > 0) {
                    pair<EdgeId, int> checking_pair = edge_stack.top();
                    edge_stack.pop();
                    if (graph_.IncomingEdgeCount(graph_.EdgeEnd(checking_pair.first)) == 1) {
                        VertexId v = graph_.EdgeEnd(checking_pair.first);
                        if (graph_.OutgoingEdgeCount(v)) {
                            for (auto I = graph_.out_begin(v), E = graph_.out_end(v); I != E; ++I) {
                                EdgeId Cur_edge = *I;
                                InTipMap.insert(
                                        std::make_pair(Cur_edge,
                                                       std::make_pair(edge,
                                                                      graph_.length(checking_pair.first) +
                                                                      checking_pair.second)));
                                edge_stack.push(
                                        std::make_pair(Cur_edge,
                                                       graph_.length(checking_pair.first) + checking_pair.second));

                            }
                        }
                    }
                }
            }

            if (graph_.OutgoingEdgeCount(graph_.EdgeEnd(edge)) == 0) {
                OutTipMap.insert(std::make_pair(edge, std::make_pair(edge, 0)));
                edge_stack.push(std::make_pair(edge, 0));
                while (edge_stack.size() > 0) {
                    std::pair<EdgeId, int> checking_pair = edge_stack.top();
                    edge_stack.pop();
                    if (graph_.OutgoingEdgeCount(graph_.EdgeStart(checking_pair.first)) == 1) {
                        if (graph_.IncomingEdgeCount(graph_.EdgeStart(checking_pair.first))) {
                            for (EdgeId e : graph_.IncomingEdges(graph_.EdgeStart(checking_pair.first))) {
                                OutTipMap.insert(std::make_pair(e,
                                                                std::make_pair(edge,
                                                                               graph_.length(e) +
                                                                               checking_pair.second)));
                                edge_stack.push(std::make_pair(e,
                                                               graph_.length(e) + checking_pair.second));
                            }
                        }
                    }

                }
            }
            ++iterator;
        }
    }

    template<class Streams>
    void MapReads(omnigraph::de::PairedInfoIndexT<Graph> &paired_index, Streams &streams,
                  const std::unordered_map<EdgeId, pair<EdgeId, int> > &OutTipMap,
                  const std::unordered_map<EdgeId, pair<EdgeId, int> > &InTipMap) const {
        INFO("Processing paired reads (takes a while)");

        size_t nthreads = streams.size();
        omnigraph::de::PairedInfoBuffersT<Graph> buffer_pi(graph_, nthreads);

        size_t counter = 0;
#       pragma omp parallel for num_threads(nthreads) reduction(+ : counter)
        for (size_t i = 0; i < nthreads; ++i) {
            typename Streams::ReadT r;
            auto &stream = streams[i];
            stream.reset();

            while (!stream.eof()) {
                stream >> r;
                ++counter;
                ProcessPairedRead(buffer_pi[i], r, OutTipMap, InTipMap);
            }
        }

        INFO("Used " << counter << " paired reads");

        INFO("Merging paired indices");
        for (auto &index: buffer_pi) {
            paired_index.Merge(index);
            index.clear();
        }
    }

public:

    GapCloserPairedIndexFiller(const Graph &graph, const SequenceMapper<Graph> &mapper)
            : graph_(graph), mapper_(mapper) { }

    /**
     * Method reads paired data from stream, maps it to genome and stores it in this PairInfoIndex.
     */
    template<class Streams>
    void FillIndex(omnigraph::de::PairedInfoIndexT<Graph> &paired_index, Streams &streams) {
        std::unordered_map<EdgeId, pair<EdgeId, int> > OutTipMap, InTipMap;

        INFO("Preparing shift maps");
        PrepareShiftMaps(OutTipMap, InTipMap);

        MapReads(paired_index, streams, OutTipMap, InTipMap);
    }

};

class GapCloser {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    Graph &g_;
    int k_;
    omnigraph::de::PairedInfoIndexT<Graph> &tips_paired_idx_;
    const size_t min_intersection_;
    const size_t hamming_dist_bound_;
    const omnigraph::de::DEWeight weight_threshold_;

    std::vector<size_t> DiffPos(const Sequence &s1, const Sequence &s2) const {
        VERIFY(s1.size() == s2.size());
        std::vector<size_t> answer;
        for (size_t i = 0; i < s1.size(); ++i)
            if (s1[i] != s2[i])
                answer.push_back(i);
        return answer;
    }

    size_t HammingDistance(const Sequence &s1, const Sequence &s2) const {
        VERIFY(s1.size() == s2.size());
        size_t dist = 0;
        for (size_t i = 0; i < s1.size(); ++i)
            if (s1[i] != s2[i])
                dist++;
        return dist;
    }

    //  size_t HammingDistance(const Sequence& s1, const Sequence& s2) const {
    //    return DiffPos(s1, s2).size();
    //  }

    vector<size_t> PosThatCanCorrect(size_t overlap_length/*in nucls*/,
                                     const vector<size_t> &mismatch_pos, size_t edge_length/*in nucls*/,
                                     bool left_edge) const {
        TRACE("Try correct left edge " << left_edge);
        TRACE("Overlap length " << overlap_length);
        TRACE("Edge length " << edge_length);
        TRACE("Mismatches " << mismatch_pos);

        vector<size_t> answer;
        for (size_t i = 0; i < mismatch_pos.size(); ++i) {
            size_t relative_mm_pos =
                    left_edge ?
                    mismatch_pos[i] :
                    overlap_length - 1 - mismatch_pos[i];
            if (overlap_length - relative_mm_pos + g_.k() < edge_length)
                //can correct mismatch
                answer.push_back(mismatch_pos[i]);
        }
        TRACE("Can correct mismatches: " << answer);
        return answer;
    }

    //todo write easier
    bool CanCorrectLeft(EdgeId e, int overlap, const vector<size_t> &mismatch_pos) const {
        return PosThatCanCorrect(overlap, mismatch_pos, g_.length(e) + g_.k(), true).size() == mismatch_pos.size();
    }

    //todo write easier
    bool CanCorrectRight(EdgeId e, int overlap,
                         const vector<size_t> &mismatch_pos) const {
        return PosThatCanCorrect(overlap, mismatch_pos, g_.length(e) + g_.k(), false).size() == mismatch_pos.size();
    }

    bool MatchesEnd(const Sequence &long_seq, const Sequence &short_seq, bool from_begin) const {
        return from_begin ? long_seq.First(short_seq.size()) == short_seq
                          : long_seq.Last(short_seq.size()) == short_seq;
    }

    void CorrectLeft(EdgeId first, EdgeId second, int overlap, const vector<size_t> &diff_pos) {
        DEBUG("Can correct first with sequence from second.");
        Sequence new_sequence = g_.EdgeNucls(first).Subseq(g_.length(first) - overlap + diff_pos.front(),
                                                           g_.length(first) + k_ - overlap)
                                + g_.EdgeNucls(second).First(k_);
        DEBUG("Checking new k+1-mers.");
        DEBUG("Check ok.");
        DEBUG("Splitting first edge.");
        pair<EdgeId, EdgeId> split_res = g_.SplitEdge(first, g_.length(first) - overlap + diff_pos.front());
        first = split_res.first;
        tips_paired_idx_.Remove(split_res.second);
        DEBUG("Adding new edge.");
        VERIFY(MatchesEnd(new_sequence, g_.VertexNucls(g_.EdgeEnd(first)), true));
        VERIFY(MatchesEnd(new_sequence, g_.VertexNucls(g_.EdgeStart(second)), false));
        g_.AddEdge(g_.EdgeEnd(first), g_.EdgeStart(second),
                new_sequence);
    }

    void CorrectRight(EdgeId first, EdgeId second, int overlap, const vector<size_t> &diff_pos) {
        DEBUG("Can correct second with sequence from first.");
        Sequence new_sequence =
                g_.EdgeNucls(first).Last(k_) + g_.EdgeNucls(second).Subseq(overlap, diff_pos.back() + 1 + k_);
        DEBUG("Checking new k+1-mers.");
        DEBUG("Check ok.");
        DEBUG("Splitting second edge.");
        pair<EdgeId, EdgeId> split_res = g_.SplitEdge(second, diff_pos.back() + 1);
        second = split_res.second;
        tips_paired_idx_.Remove(split_res.first);
        DEBUG("Adding new edge.");
        VERIFY(MatchesEnd(new_sequence, g_.VertexNucls(g_.EdgeEnd(first)), true));
        VERIFY(MatchesEnd(new_sequence, g_.VertexNucls(g_.EdgeStart(second)), false));

        g_.AddEdge(g_.EdgeEnd(first), g_.EdgeStart(second),
                new_sequence);
    }

    bool HandlePositiveHammingDistanceCase(EdgeId first, EdgeId second, int overlap) {
        DEBUG("Match was imperfect. Trying to correct one of the tips");
        vector<size_t> diff_pos = DiffPos(g_.EdgeNucls(first).Last(overlap),
                                          g_.EdgeNucls(second).First(overlap));
        if (CanCorrectLeft(first, overlap, diff_pos)) {
            CorrectLeft(first, second, overlap, diff_pos);
            return true;
        } else if (CanCorrectRight(second, overlap, diff_pos)) {
            CorrectRight(first, second, overlap, diff_pos);
            return true;
        } else {
            DEBUG("Can't correct tips due to the graph structure");
            return false;
        }
    }

    bool HandleSimpleCase(EdgeId first, EdgeId second, int overlap) {
        DEBUG("Match was perfect. No correction needed");
        DEBUG("Overlap " << overlap);
        //strange info guard
        VERIFY(overlap <= k_);
        if (overlap == k_) {
            DEBUG("Tried to close zero gap");
            return false;
        }
        //old code
        Sequence edge_sequence = g_.EdgeNucls(first).Last(k_)
                                 + g_.EdgeNucls(second).Subseq(overlap, k_);
        DEBUG("Gap filled: Gap size = " << k_ - overlap << "  Result seq "
              << edge_sequence.str());
        g_.AddEdge(g_.EdgeEnd(first), g_.EdgeStart(second), edge_sequence);
        return true;
    }

    bool ProcessPair(EdgeId first, EdgeId second) {
        TRACE("Processing edges " << g_.str(first) << " and " << g_.str(second));
        TRACE("first " << g_.EdgeNucls(first) << " second " << g_.EdgeNucls(second));

        if (cfg::get().avoid_rc_connections &&
            (first == g_.conjugate(second) || first == second)) {
            DEBUG("Trying to join conjugate edges " << g_.int_id(first));
            return false;
        }

        Sequence seq1 = g_.EdgeNucls(first);
        Sequence seq2 = g_.EdgeNucls(second);
        TRACE("Checking possible gaps from 1 to " << k_ - min_intersection_);
        for (int gap = 1; gap <= k_ - (int) min_intersection_; ++gap) {
            int overlap = k_ - gap;
            size_t hamming_distance = HammingDistance(g_.EdgeNucls(first).Last(overlap),
                                                      g_.EdgeNucls(second).First(overlap));
            if (hamming_distance <= hamming_dist_bound_) {
                DEBUG("For edges " << g_.str(first) << " and " << g_.str(second)
                      << ". For gap value " << gap << " (overlap " << overlap << "bp) hamming distance was " <<
                      hamming_distance);
                //        DEBUG("Sequences of distance " << tip_distance << " :"
                //                << seq1.Subseq(seq1.size() - k).str() << "  "
                //                << seq2.Subseq(0, k).str());

                if (hamming_distance > 0) {
                    return HandlePositiveHammingDistanceCase(first, second, overlap);
                } else {
                    return HandleSimpleCase(first, second, overlap);
                }
            }
        }
        return false;
    }

public:
    //TODO extract methods
    void CloseShortGaps() {
        INFO("Closing short gaps");
        size_t gaps_filled = 0;
        size_t gaps_checked = 0;
        for (auto edge = g_.SmartEdgeBegin(); !edge.IsEnd(); ++edge) {
            EdgeId first_edge = *edge;
            for (auto i : tips_paired_idx_.Get(first_edge)) {
                EdgeId second_edge = i.first;
                if (first_edge == second_edge)
                    continue;

                if (!g_.IsDeadEnd(g_.EdgeEnd(first_edge)) || !g_.IsDeadStart(g_.EdgeStart(second_edge))) {
                    // WARN("Topologically wrong tips");
                    continue;
                }

                bool closed = false;
                for (auto point : i.second) {
                    if (math::ls(point.weight, weight_threshold_))
                        continue;

                    ++gaps_checked;
                    closed = ProcessPair(first_edge, second_edge);
                    if (closed) {
                        ++gaps_filled;
                        break;
                    }
                }
                if (closed)
                    break;
            } // second edge
        } // first edge

        INFO("Closing short gaps complete: filled " << gaps_filled
             << " gaps after checking " << gaps_checked
             << " candidates");
        omnigraph::CompressAllVertices<Graph>(g_);
    }

    GapCloser(Graph &g, omnigraph::de::PairedInfoIndexT<Graph> &tips_paired_idx,
              size_t min_intersection, double weight_threshold,
              size_t hamming_dist_bound = 0 /*min_intersection_ / 5*/)
            : g_(g),
              k_((int) g_.k()),
              tips_paired_idx_(tips_paired_idx),
              min_intersection_(min_intersection),
              hamming_dist_bound_(hamming_dist_bound),
              weight_threshold_(weight_threshold)  {
        VERIFY(min_intersection_ < g_.k());
        DEBUG("weight_threshold=" << weight_threshold_);
        DEBUG("min_intersect=" << min_intersection_);
        DEBUG("paired_index size=" << tips_paired_idx_.size());
    }

private:
    DECL_LOGGER("GapCloser");
};

template<class Streams>
void CloseGaps(conj_graph_pack &gp, Streams &streams) {
    auto mapper = MapperInstance(gp);
    GapCloserPairedIndexFiller gcpif(gp.g, *mapper);
    PairedIndexT tips_paired_idx(gp.g);
    gcpif.FillIndex(tips_paired_idx, streams);
    GapCloser gap_closer(gp.g, tips_paired_idx,
                         cfg::get().gc.minimal_intersection, cfg::get().gc.weight_threshold);
    gap_closer.CloseShortGaps();
}

void GapClosing::run(conj_graph_pack &gp, const char *) {
    visualization::graph_labeler::DefaultLabeler<Graph> labeler(gp.g, gp.edge_pos);
    stats::detail_info_printer printer(gp, labeler, cfg::get().output_dir);
    printer(config::info_printer_pos::before_first_gap_closer);

    bool pe_exist = false;
    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
        if (cfg::get().ds.reads[i].type() == io::LibraryType::PairedEnd) {
            pe_exist = true;
            break;
        }
    }
    if (!pe_exist) {
        INFO("No paired-end libraries exist, skipping gap closer");
        return;
    }
    gp.EnsureIndex();

    auto& dataset = cfg::get_writable().ds;
    for (size_t i = 0; i < dataset.reads.lib_count(); ++i) {
        if (dataset.reads[i].type() == io::LibraryType::PairedEnd) {
            auto streams = paired_binary_readers(dataset.reads[i], false, 0, false);
            CloseGaps(gp, streams);
        }
    }
}

}
