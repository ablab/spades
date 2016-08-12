//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/graph_core/graph.hpp"
#include "assembly_graph/graph_alignment/sequence_mapper.hpp"
#include "ConsensusCore/Poa/PoaConfig.hpp"
#include "ConsensusCore/Poa/PoaConsensus.hpp"
#include "gap_closing.hpp"

#include <algorithm>
#include <fstream>

namespace debruijn_graph {
namespace gap_closing {

class GapStorage {
public:
    typedef vector<GapDescription> GapInfos;
    typedef typename GapInfos::const_iterator gap_info_it;
private:

    const Graph& g_;
    const size_t min_gap_quantity_;
    const size_t long_seq_limit_;
    const size_t max_flanking_region_length_;

    map<EdgeId, GapInfos> inner_index_;
    vector<EdgeId> index_;
    set<pair<EdgeId, EdgeId>> transitively_ignored_pairs_;
    set<pair<EdgeId, EdgeId>> symmetrically_ignored_pairs_;

    DECL_LOGGER("GapStorage");

    void HiddenAddGap(const GapDescription& p) {
        inner_index_[p.start].push_back(p);
    }

    bool CheckGap(const GapDescription& gap) const {
        return gap.edge_gap_start_position + max_flanking_region_length_ > g_.length(gap.start)
               && gap.edge_gap_end_position < max_flanking_region_length_;
    }

    //all gaps guaranteed to correspond to a single edge pair
    GapInfos PadEdgePairGaps(const gap_info_it& start, const gap_info_it& end) const {
        size_t start_min = std::numeric_limits<size_t>::max();
        size_t end_max = 0;
        size_t long_seqs = 0;
        size_t short_seqs = 0;
        for (auto it = start; it != end; ++it) {
            const auto& gap = *it;
            if (CheckGap(gap)) {
                if (gap.gap_seq.size() > long_seq_limit_)
                    long_seqs++;
                else
                    short_seqs++;

                start_min = std::min(start_min, gap.edge_gap_start_position);
                end_max = std::max(end_max, gap.edge_gap_end_position);
            } else {
                DEBUG("ignoring alingment to the middle of edge");
            }
        }

        const bool exclude_long_seqs = (short_seqs >= min_gap_quantity_ && short_seqs > long_seqs);

        GapInfos answer;
        for (auto it = start; it != end; ++it) {
            const auto& gap = *it;
            if (!CheckGap(gap))
                continue;

            if (exclude_long_seqs && gap.gap_seq.size() > long_seq_limit_)
                continue;

            string s = g_.EdgeNucls(gap.start).Subseq(start_min, gap.edge_gap_start_position).str();
            s += gap.gap_seq.str();
            s += g_.EdgeNucls(gap.end).Subseq(gap.edge_gap_end_position, end_max).str();
            answer.push_back(GapDescription(gap.start, gap.end, Sequence(s), start_min, end_max));
        }
        return answer;
    }

    GapInfos PadEdgeGaps(const GapInfos& edge_gaps) const {
        GapInfos answer;

        for (const auto& edge_pair_gaps: EdgePairGaps(edge_gaps)) {
            push_back_all(answer, PadEdgePairGaps(edge_pair_gaps.first, edge_pair_gaps.second));
        }
        return answer;
    }

    size_t FillIndex() {
        VERIFY(index_.empty());
        index_.reserve(inner_index_.size());
        set<EdgeId> tmp;
        for (const auto& kv : inner_index_) {
            index_.push_back(kv.first);
        }
        return index_.size();
    }


public:

    GapStorage(const Graph& g, size_t min_gap_quantity, size_t long_seq_limit,
               size_t max_flanking_region_length = 500)
            : g_(g),
              min_gap_quantity_(min_gap_quantity),
              long_seq_limit_(long_seq_limit),
              max_flanking_region_length_(max_flanking_region_length) {
    }

    const map<EdgeId, GapInfos>& inner_index() const {
        return inner_index_;
    };

    size_t min_gap_quantity() const {
        return min_gap_quantity_;
    }

    EdgeId operator[](size_t i) const {
        return index_.at(i);
    }

    size_t size() const {
        return index_.size();
    }

    bool IsIgnored(const pair<EdgeId, EdgeId>& p) const {
        return transitively_ignored_pairs_.count(p) || symmetrically_ignored_pairs_.count(p);
    }

    void AddGap(const GapDescription& p, bool add_rc = false) {
        HiddenAddGap(p);
        if (add_rc) {
            HiddenAddGap(p.conjugate(g_));
        }
    }

    void AddStorage(const GapStorage& to_add) {
        const auto& idx = to_add.inner_index_;
        for (auto iter = idx.begin(); iter != idx.end(); ++iter)
            inner_index_[iter->first].insert(inner_index_[iter->first].end(), iter->second.begin(), iter->second.end());
    }

    void clear() {
        GapStorage empty(g_, min_gap_quantity_, long_seq_limit_);
        std::swap(inner_index_, empty.inner_index_);
        std::swap(index_, empty.index_);
        std::swap(transitively_ignored_pairs_, empty.transitively_ignored_pairs_);
        std::swap(symmetrically_ignored_pairs_, empty.symmetrically_ignored_pairs_);
    }

    void DumpToFile(const string filename) const {
        ofstream filestr(filename);
        for (const auto& e_gaps : inner_index_) {
            EdgeId e = e_gaps.first;
            auto gaps = e_gaps.second;
            DEBUG(g_.int_id(e) << " " << gaps.size());
            filestr << g_.int_id(e) << " " << gaps.size() << endl;
            sort(gaps.begin(), gaps.end());
            for (const auto& gap : gaps) {
                filestr << gap.str(g_);
            }
            filestr << endl;
        }
    }

//    void LoadFromFile(const string s) {
//        FILE* file = fopen((s).c_str(), "r");
//        int res;
//        char ss[5000];
//        map<int, EdgeId> tmp_map;
//        for (auto iter = g.ConstEdgeBegin(); !iter.IsEnd(); ++iter) {
//            tmp_map[g.int_id(*iter)] = *iter;
//        }
//        while (!feof(file)) {
//            int first_id, second_id, first_ind, second_ind;
//            int size;
//            res = fscanf(file, "%d %d\n", &first_id, &size);
//            VERIFY(res == 2);
//            for (int i = 0; i < size; i++) {
//                res = fscanf(file, "%d %d\n", &first_id, &first_ind);
//                VERIFY(res == 2);
//                res = fscanf(file, "%d %d\n", &second_id, &second_ind);
//                VERIFY(res == 2);
//                res = fscanf(file, "%s\n", ss);
//                VERIFY(res == 1);
//                GapDescription<Graph> gap(tmp_map[first_id], tmp_map[second_id], Sequence(ss), first_ind, second_ind);
//                this->AddGap(gap);
//            }
//        }
//    }

    //edge_gaps must be sorted
    vector<pair<gap_info_it, gap_info_it>> EdgePairGaps(const GapInfos& edge_gaps) const {
        vector<pair<gap_info_it, gap_info_it>> answer;
        auto ep_start = edge_gaps.begin();
        for (auto it = ep_start; it != edge_gaps.end(); ++it) {
            if (it->end != ep_start->end) {
                answer.push_back({ep_start, it});
                ep_start = it;
            }
        }
        answer.push_back({ep_start, edge_gaps.end()});
        return answer;
    };


    void PrepareGapsForClosure() {
        for (auto& e_gaps : inner_index_) {
            DEBUG("Padding gaps for first edge " << g_.str(e_gaps.first));
            auto& gaps = e_gaps.second;
            sort(gaps.begin(), gaps.end());
            gaps = PadEdgeGaps(gaps);
        }

        FillIndex();

        set<pair<EdgeId, EdgeId>> nonempty_pairs;
        for (auto& e_gaps : inner_index_) {
            for (const auto& edge_pair_gaps: EdgePairGaps(e_gaps.second)) {
                if (edge_pair_gaps.second - edge_pair_gaps.first >= min_gap_quantity_) {
                    nonempty_pairs.insert({edge_pair_gaps.first->start, edge_pair_gaps.first->end});
                }
            }
        }

        set<pair<EdgeId, EdgeId>> used_rc_pairs;
        for (const auto& edge_pair : nonempty_pairs) {
            if (used_rc_pairs.count(edge_pair)) {
                DEBUG("skipping pair " << g_.int_id(edge_pair.first) << "," << g_.int_id(edge_pair.second));
                symmetrically_ignored_pairs_.insert(edge_pair);
            } else {
                DEBUG("Using pair " << g_.int_id(edge_pair.first) << "," << g_.int_id(edge_pair.second));
            }

            for (EdgeId e : index_) {
                if (nonempty_pairs.count(make_pair(edge_pair.first, e))
                    && nonempty_pairs.count(make_pair(e, edge_pair.second))) {
                    DEBUG("pair " << g_.int_id(edge_pair.first) << "," << g_.int_id(edge_pair.second)
                                  << " is ignored because of edge between " << g_.int_id(e));
                    transitively_ignored_pairs_.insert(edge_pair);
                }
            }
            used_rc_pairs.insert(make_pair(g_.conjugate(edge_pair.second), g_.conjugate(edge_pair.first)));
        }
    }
};

inline string PoaConsensus(const vector<string>& gap_seqs) {
    const ConsensusCore::PoaConsensus* pc = ConsensusCore::PoaConsensus::FindConsensus(
            gap_seqs,
            ConsensusCore::PoaConfig::GLOBAL_ALIGNMENT);
    return pc->Sequence();
}

inline string TrivialConsenus(const vector<string>& gap_seqs, size_t max_length) {
    VERIFY(!gap_seqs.empty());
    return gap_seqs.front().length() < max_length ? gap_seqs.front() : "";
}

class HybridGapCloser {
public:
    typedef std::function<string (const vector<string>&)> ConsensusF;
private:
    typedef runtime_k::RtSeq Kmer;
    typedef typename GapStorage::gap_info_it gap_info_it;

    DECL_LOGGER("HybridGapCloser");

    Graph& g_;
    const GapStorage& storage_;
    const size_t min_weight_;
    ConsensusF consensus_;

    const GapDescription INVALID_GAP;

    string PrintLengths(const vector<string>& gap_seqs) const {
        stringstream ss;
        for (const auto& gap_v : gap_seqs)
            ss << gap_v.length() << " ";
        return ss.str();
    }

    GapDescription ConstructConsensus(EdgeId start,
                                      EdgeId end,
                                      size_t edge_gap_start_position,
                                      size_t edge_gap_end_position,
                                      const vector<string>& gap_variants) const {
        if (gap_variants.size() > 1) {
            DEBUG(gap_variants.size() << " gap closing variants, lengths: " << PrintLengths(gap_variants));
        }
        string s = consensus_(gap_variants);
        if (!s.empty()) {
            DEBUG("consenus for " << g_.int_id(start)
                                  << " and " << g_.int_id(end)
                                  << "found: " << s);
            return GapDescription(start, end,
                                  Sequence(s),
                                  edge_gap_start_position, edge_gap_end_position);
        } else {
            INFO("Skipping gap of size " << gap_variants.front().length() << " multiplicity " << gap_variants.size());
        }

        return INVALID_GAP;
    }

    GapDescription ConstructConsensus(gap_info_it start_it, gap_info_it end_it) const {
        size_t cur_len = end_it - start_it;

        if (cur_len < min_weight_ || storage_.IsIgnored(make_pair(start_it->start, start_it->end)))
            return INVALID_GAP;

        vector<string> gap_variants;
        std::transform(start_it, end_it, std::back_inserter(gap_variants), [](const GapDescription& gap) {
            return gap.gap_seq.str();
        });

        //for (auto it = start_it; it != end_it; ++it) {
        //    VERIFY(it->start == start_it->start);
        //    VERIFY(it->end == start_it->end);
        //    VERIFY(it->edge_gap_start_position == start_it->edge_gap_start_position);
        //    VERIFY(it->edge_gap_end_position == start_it->edge_gap_end_position);
        //}

        return ConstructConsensus(start_it->start, start_it->end,
                                  start_it->edge_gap_start_position,
                                  start_it->edge_gap_end_position,
                                  gap_variants);
    }

    GapDescription ConstructConsensus(EdgeId e) const {
        vector<GapDescription> closures;
        for (const auto& edge_pair_gaps : storage_.EdgePairGaps(get(storage_.inner_index(), e))) {
            auto consensus = ConstructConsensus(edge_pair_gaps.first, edge_pair_gaps.second);
            if (consensus != INVALID_GAP) {
                closures.push_back(consensus);
            }
        }
        if (closures.size() == 1)
            return closures.front();

        if (closures.size() > 1)
            DEBUG("non-unique gap!!");
        return INVALID_GAP;
    }

    vector<GapDescription> ConstructConsensus(size_t nthreads) const {
        vector<vector<GapDescription>> closures_by_thread(nthreads);

        # pragma omp parallel for num_threads(nthreads)
        for (size_t i = 0; i < storage_.size(); i++) {
            EdgeId e = storage_[i];
            size_t thread_num = omp_get_thread_num();
            DEBUG("constructing consenus for first edge " << g_.int_id(e) << " in thread " << thread_num);
            GapDescription gap = ConstructConsensus(e);
            if (gap != INVALID_GAP) {
                closures_by_thread[thread_num].push_back(gap);
            }
        }

        vector<GapDescription> closures;
        for (auto& new_per_thread : closures_by_thread) {
            std::copy(new_per_thread.begin(), new_per_thread.end(), std::back_inserter(closures));
            new_per_thread.clear();
        }
        return closures;
    }

    map<EdgeId, EdgeId> CloseGapsInGraph(const vector<GapDescription>& closures) {
        map<EdgeId, EdgeId> replacement;
        size_t closed_gaps = 0;
        GapJoiner gap_joiner(g_);
        for (auto new_edge_info : closures) {
            EdgeId first = new_edge_info.start;
            EdgeId second = new_edge_info.end;
            if (replacement.count(first) || replacement.count(second)) {
                DEBUG("sorry, gap chains are not supported yet");
                continue;
            }

            EdgeId new_edge = gap_joiner(new_edge_info, /*compress*/true);

            TRACE("New edge " << g_.str(new_edge) << " coverage " << g_.coverage(new_edge));
            closed_gaps++;

            replacement[first] = new_edge;
            replacement[second] = new_edge;
            replacement[g_.conjugate(first)] = g_.conjugate(new_edge);
            replacement[g_.conjugate(second)] = g_.conjugate(new_edge);
        }
        INFO("Closed " << closed_gaps << " gaps");
        //TODO: chains of gaps!
        return replacement;
    }

public:
    HybridGapCloser(Graph& g, const GapStorage& storage,
                    size_t min_weight, ConsensusF consensus)
            : g_(g), storage_(storage),
              min_weight_(min_weight),
              consensus_(consensus) {
    }

    map<EdgeId, EdgeId> operator()(size_t nthreads) {
        return CloseGapsInGraph(ConstructConsensus(nthreads));
    };

};

}
}
