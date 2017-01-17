//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/core/graph.hpp"
#include "modules/alignment/sequence_mapper.hpp"
#include "ConsensusCore/Poa/PoaConfig.hpp"
#include "ConsensusCore/Poa/PoaConsensus.hpp"
#include "gap_closing.hpp"

#include <algorithm>
#include <fstream>

namespace debruijn_graph {
namespace gap_closing {
typedef vector<GapDescription> GapInfos;

typedef pair<EdgeId, EdgeId> EdgePair;
inline EdgePair Conjugate(const Graph& g, EdgePair ep) {
    return EdgePair(g.conjugate(ep.second), g.conjugate(ep.first));
}

inline bool IsCanonical(const Graph& g, const EdgePair& ep) {
    return ep <= Conjugate(g, ep);
}

inline bool IsCanonical(const Graph& g, EdgeId a, EdgeId b) {
    return IsCanonical(g, EdgePair(a,b));
}

inline bool IsCanonical(const Graph& g, EdgeId e) {
    return e <= g.conjugate(e);
}

inline EdgePair GetCanonical(const Graph& g, const EdgePair& ep) {
    return IsCanonical(g, ep) ? ep : Conjugate(g, ep);
}

class GapStorage {
public:
    typedef typename GapInfos::const_iterator gap_info_it;
    typedef std::pair<gap_info_it, gap_info_it> info_it_pair;
private:
    typedef std::function<bool (gap_info_it, gap_info_it)> CandidatesPred;
    typedef std::function<bool (const EdgePair&)> EdgePairPred;
    typedef std::function<bool (const GapDescription&)> DescriptionPred;
    typedef std::set<EdgePair> ConnectionSet;

    const Graph& g_;

    map<EdgeId, GapInfos> inner_index_;
    vector<EdgeId> index_;

    DECL_LOGGER("GapStorage");

    void HiddenAddGap(const GapDescription& p) {
        inner_index_[p.start].push_back(p);
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

    typename std::vector<GapDescription>::iterator
    const_iterator_cast(std::vector<GapDescription> &v,
                        typename std::vector<GapDescription>::const_iterator iter) const {
        return v.begin() + (iter - v.cbegin());
    }

    //Function should return true if corresponding part of the index should be removed
    void FilterByCandidates(const CandidatesPred &filter_f) {
        for (auto it = inner_index_.begin(); it != inner_index_.end(); ) {
            vector<GapDescription>& gaps = it->second;
            auto ep_ranges = EdgePairGaps(gaps);

            auto copy_dest = gaps.begin();
            for (const info_it_pair& ep_gaps : ep_ranges) {
                if (filter_f(ep_gaps.first, ep_gaps.second)) {
                    DEBUG("Erasing candidates between " << g_.int_id(ep_gaps.first->start) << " and "
                                                        << g_.int_id(ep_gaps.first->end));
                } else {
                    if (copy_dest == const_iterator_cast(gaps, ep_gaps.first)) {
                        copy_dest = const_iterator_cast(gaps, ep_gaps.second);
                    } else {
                        copy_dest = std::move(ep_gaps.first, ep_gaps.second, copy_dest);
                    }
                }
            }
            if (copy_dest == gaps.begin()) {
                inner_index_.erase(it++);
            } else {
                gaps.erase(copy_dest, gaps.end());
                ++it;
            }
        }
    }

    void FilterByEdgePair(const EdgePairPred &filter_f) {
        FilterByCandidates([=](gap_info_it info_start, gap_info_it /*info_end*/) {
            return filter_f(EdgePair(info_start->start, info_start->end));
        });
    }

    void FilterByDescription(const DescriptionPred &filter_f) {
        for (auto it = inner_index_.begin(); it != inner_index_.end(); ) {
            vector<GapDescription>& gaps = it->second;
            auto res_it = std::remove_if(gaps.begin(), gaps.end(), filter_f);
            if (res_it == gaps.begin()) {
                inner_index_.erase(it++);
            } else {
                gaps.erase(res_it, gaps.end());
                ++it;
            }
        }
    }

    vector<EdgeId> SecondEdges(const GapInfos& edge_gaps) const {
        vector<EdgeId> jump_edges;
        for (auto it_pair : EdgePairGaps(edge_gaps)) {
            jump_edges.push_back(it_pair.first->end);
        }
        return jump_edges;
    };

    ConnectionSet GetAllConnections() const {
        ConnectionSet answer;
        for (const auto& e_gaps : inner_index_) {
            EdgeId e1 = e_gaps.first;
            for (EdgeId e2: SecondEdges(e_gaps.second)) {
                EdgePair ep(e1, e2);
                answer.insert(ep);
                answer.insert(Conjugate(g_, ep));
            }
        }
        return answer;
    };

    //outputs set of transitively-redundant CANONICAL connections
    ConnectionSet DetectTransitive() const {
        auto all_connections = GetAllConnections();
        ConnectionSet answer;
        for (auto it = all_connections.begin(), end_it = all_connections.end(); it != end_it; ) {
            EdgeId left = it->first;
            vector<EdgeId> right_options;
            auto inner_it = it;
            for (; inner_it != end_it && inner_it->first == left; ++inner_it) {
                right_options.push_back(inner_it->second);
            }

            for (size_t i = 0; i < right_options.size(); ++i) {
                for (size_t j = 0; j < right_options.size(); ++j) {
                    if (i == j)
                        continue;
                    if (all_connections.count(EdgePair(right_options[i], right_options[j]))) {
                        //TODO should we add sanity checks that other edges of the triangle are not there?
                        answer.insert(GetCanonical(g_, EdgePair(left, right_options[j])));
                        DEBUG("pair " << g_.int_id(left) << "," << g_.int_id(right_options[j])
                                      << " is ignored because of edge between "
                                      << g_.int_id(right_options[i]));
                    }
                }
            }
            it = inner_it;
        }
        return answer;
    }

    std::set<EdgeId> AmbiguouslyExtending() const {
        std::set<EdgeId> answer;
        std::set<EdgeId> left_edges;
        for (const auto& e_gaps : inner_index_) {
            EdgeId e1 = e_gaps.first;
            for (EdgeId e2: SecondEdges(e_gaps.second)) {
                if (!left_edges.insert(e1).second) {
                    answer.insert(e1);
                }
                if (!left_edges.insert(g_.conjugate(e2)).second) {
                    answer.insert(g_.conjugate(e2));
                }
            }
        }
        return answer;
    }

    void FilterIndex(size_t min_weight, size_t max_flank) {
        DEBUG("Filtering by maximal allowed flanking length " << max_flank);
        FilterByDescription([=](const GapDescription &gap) {
            return gap.edge_gap_start_position + max_flank < g_.length(gap.start)
                   || gap.edge_gap_end_position > max_flank;
        });

        DEBUG("Filtering by weight " << min_weight);
        FilterByCandidates([=](gap_info_it info_start, gap_info_it info_end) {
            auto cnt = std::distance(info_start, info_end);
            VERIFY(cnt > 0);
            return size_t(cnt) < min_weight;
        });


        DEBUG("Filtering transitive gaps");
        ConnectionSet transitive_ignore = DetectTransitive();

        FilterByEdgePair([&](const EdgePair &ep) {
            VERIFY(IsCanonical(g_, ep));
            return transitive_ignore.count(ep);
        });

        DEBUG("Filtering ambiguous situations");
        std::set<EdgeId> ambiguously_extending = AmbiguouslyExtending();
        FilterByEdgePair([&](const EdgePair &ep) {
            return ambiguously_extending.count(ep.first) ||
                    ambiguously_extending.count(g_.conjugate(ep.second));
        });
    }

public:

    GapStorage(const Graph& g)
            : g_(g) {
    }

    const map<EdgeId, GapInfos>& inner_index() const {
        return inner_index_;
    };

    EdgeId operator[](size_t i) const {
        return index_.at(i);
    }

    size_t size() const {
        return index_.size();
    }

    void AddGap(const GapDescription& p) {
        if (IsCanonical(g_, p.start, p.end)) {
            HiddenAddGap(p);
        } else {
            HiddenAddGap(p.conjugate(g_));
        }
    }

    void AddStorage(const GapStorage& to_add) {
        const auto& idx = to_add.inner_index_;
        for (auto iter = idx.begin(); iter != idx.end(); ++iter)
            inner_index_[iter->first].insert(inner_index_[iter->first].end(), iter->second.begin(), iter->second.end());
    }

    void clear() {
        GapStorage empty(g_);
        std::swap(inner_index_, empty.inner_index_);
        std::swap(index_, empty.index_);
    }

    void DumpToFile(const string filename) const {
        ofstream filestr(filename);
        for (const auto& e_gaps : inner_index_) {
            EdgeId e = e_gaps.first;
            auto gaps = e_gaps.second;
            DEBUG(g_.int_id(e) << " " << gaps.size());
            filestr << g_.int_id(e) << " " << gaps.size() << endl;
            std::sort(gaps.begin(), gaps.end());
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
    vector<info_it_pair> EdgePairGaps(const GapInfos& edge_gaps) const {
        vector<info_it_pair> answer;
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

    void PrepareGapsForClosure(size_t min_weight, size_t max_flank) {
        for (auto& e_gaps : inner_index_) {
            auto& gaps = e_gaps.second;
            std::sort(gaps.begin(), gaps.end());
        }
        DEBUG("Raw extensions available for " << inner_index_.size() << " edges");

        FilterIndex(min_weight, max_flank);
        DEBUG("Filtered extensions available for " << inner_index_.size() << " edges");
        FillIndex();
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

/*Keys are actual edges of the graph, values are original edges*/
/*In general many-to-many relationship*/
class EdgeFateTracker : omnigraph::GraphActionHandler<Graph> {
    map<EdgeId, set<EdgeId>> storage_;

    void FillRelevant(EdgeId e, set<EdgeId>& relevant) const {
        auto it = storage_.find(e);
        if (it != storage_.end()) {
            //one of novel edges
            relevant.insert(it->second.begin(), it->second.end());
        } else {
            //one of original edges
            relevant.insert(e);
        }
    }

public:
    EdgeFateTracker(const Graph& g) :
            omnigraph::GraphActionHandler<Graph>(g, "EdgeFateTracker") {
    }

    void HandleAdd(EdgeId e) override {
        if (!storage_.count(e))
            storage_[e] = {};
    }

    void HandleDelete(EdgeId e) override {
        storage_.erase(e);
    }

    void HandleMerge(const vector<EdgeId>& old_edges, EdgeId new_edge) override {
        set<EdgeId> relevant_records;
        for (EdgeId e : old_edges) {
            FillRelevant(e, relevant_records);
        }
        storage_[new_edge] = relevant_records;
    }

    void HandleGlue(EdgeId /*new_edge*/, EdgeId /*edge1*/, EdgeId /*edge2*/) override {
        VERIFY(false);
    }

    void HandleSplit(EdgeId old_edge, EdgeId new_edge_1,
                             EdgeId new_edge_2) override {
        set<EdgeId> relevant_records;
        FillRelevant(old_edge, relevant_records);
        storage_[new_edge_1] = relevant_records;
        storage_[new_edge_2] = relevant_records;
    }

    map<EdgeId, EdgeId> Old2NewMapping() const {
        map<EdgeId, EdgeId> old_2_new;
        for (const auto& new_2_olds : storage_) {
            for (EdgeId e : new_2_olds.second) {
                VERIFY(!old_2_new.count(e));
                old_2_new[e] = new_2_olds.first;
            }
        }
        return old_2_new;
    }

};

class MultiGapJoiner {
    typedef map<EdgeId, pair<size_t, size_t>> SplitInfo;

    Graph& g_;
    GapJoiner inner_joiner_;

    bool CheckGapsValidity(const vector<GapDescription>& gaps) const {
        vector<GapDescription> answer;
        return std::all_of(gaps.begin(), gaps.end(), [&](const GapDescription &gap) {
            return IsCanonical(g_, gap.start, gap.end) && gap.start != gap.end && gap.start != g_.conjugate(gap.end);
        });
    }

    void Add(size_t idx, EdgeId e, size_t pos, SplitInfo& primary, SplitInfo& secondary) const {
        SplitInfo* storage = &primary;
        if (!IsCanonical(g_, e)) {
            e = g_.conjugate(e);
            pos = g_.length(e) - pos;
            storage = &secondary;
        }
        VERIFY(!storage->count(e));
        storage->insert(make_pair(e, make_pair(idx, pos)));
    }

    vector<EdgeId> EdgesNeedingSplit(const SplitInfo& left_split_info, const SplitInfo& right_split_info) const {
        vector<EdgeId> answer;
        for (EdgeId e : key_set(left_split_info))
            if (right_split_info.count(e))
                answer.push_back(e);
        return answer;
    }

    size_t ArtificialSplitPos(size_t left_split, size_t right_split) const {
        if (right_split < left_split + 2) {
            DEBUG("Artificial split impossible");
            return -1ul;
        }
        return (left_split + right_split) / 2;
    }

    bool Update(EdgeId& e, size_t& gap_pos, EdgePair split_orig_ep, EdgePair split_res, bool gap_start) const {
        EdgeId split_orig = split_orig_ep.first;
        if (e == split_orig_ep.second) {
            split_orig = split_orig_ep.second;
            split_res = Conjugate(g_, split_res);
        }
        if (e == split_orig) {
            if (gap_start) {
                e = split_res.second;
                gap_pos = gap_pos - g_.length(split_res.first);
            } else {
                e = split_res.first;
            }
            return true;
        }
        return false;
    }

    void UpdateGap(GapDescription& gap, EdgePair split_orig, EdgePair split_res) const {
        bool u1 = Update(gap.start, gap.edge_gap_start_position, split_orig, split_res, true);
        bool u2 = Update(gap.end, gap.edge_gap_end_position, split_orig, split_res, false);
        VERIFY(u1 != u2);
    }

    bool CheckInsert(EdgeId e, set<EdgeId>& used_edges) const {
        return used_edges.insert(e).second;
    }

    bool CheckInsert(const vector<EdgeId> edges, set<EdgeId>& used_edges) const {
        for (EdgeId e : edges) {
            if (!CheckInsert(e, used_edges)) {
                return false;
            }
        }
        return true;
    }

    std::set<EdgeId> RelevantEdges(const GapDescription& gap) const {
        std::set<EdgeId> answer;
        answer.insert(gap.start);
        answer.insert(g_.conjugate(gap.start));
        answer.insert(gap.end);
        answer.insert(g_.conjugate(gap.end));
        return answer;
    }

    bool CheckGaps(const vector<GapDescription>& gaps) const {
        set<EdgeId> used_edges;
        for (const auto& gap : gaps) {
            const auto relevant = RelevantEdges(gap);
            //TODO check the semantics of all_of
            if (!std::all_of(relevant.begin(), relevant.end(), [&](const EdgeId& e) {
                return used_edges.insert(e).second;
            })) {
                return false;
            }
        }
        return true;
    }

    vector<GapDescription> ArtificialSplitAndGapUpdate(vector<GapDescription> canonical_gaps) const {
        SplitInfo left_split_pos;
        SplitInfo right_split_pos;
        for (size_t i = 0; i < canonical_gaps.size(); ++i) {
            const auto& gap = canonical_gaps[i];
            DEBUG("Processing gap " << gap.str(g_));
            Add(i, gap.start, gap.edge_gap_start_position, right_split_pos, left_split_pos);
            Add(i, gap.end, gap.edge_gap_end_position, left_split_pos, right_split_pos);
        }

        set<size_t> to_ignore;

        for (EdgeId e : EdgesNeedingSplit(left_split_pos, right_split_pos)) {
            size_t artificial_split_pos = ArtificialSplitPos(left_split_pos[e].second, right_split_pos[e].second);
            if (artificial_split_pos == -1ul) {
                to_ignore.insert(left_split_pos[e].first);
                to_ignore.insert(right_split_pos[e].first);
            } else {
                DEBUG("Splitting edge " << g_.str(e) << " at pos " << artificial_split_pos);
                DEBUG("Will update gap " << canonical_gaps[left_split_pos[e].first].str(g_) << " and " << canonical_gaps[right_split_pos[e].first].str(g_));
                EdgePair ep(e, g_.conjugate(e));
                auto split_res = g_.SplitEdge(e, artificial_split_pos);
                UpdateGap(canonical_gaps[left_split_pos[e].first], ep, split_res);
                UpdateGap(canonical_gaps[right_split_pos[e].first], ep, split_res);
            }
        }

        vector<GapDescription> updated_gaps;
        updated_gaps.reserve(canonical_gaps.size());
        for (size_t i = 0; i < canonical_gaps.size(); ++i) {
            if (!to_ignore.count(i)) {
                updated_gaps.push_back(canonical_gaps[i]);
            }
        }

        VERIFY(CheckGaps(updated_gaps));
        return updated_gaps;
    };

public:
    MultiGapJoiner(Graph& g) : g_(g), inner_joiner_(g, true) {
    }

    //Resulting graph should be condensed
    void operator()(const vector<GapDescription>& gaps) {
        size_t closed_gaps = 0;
        VERIFY_MSG(CheckGapsValidity(gaps), "Gap check failed");
        for (const auto& gap : ArtificialSplitAndGapUpdate(gaps)) {
            inner_joiner_(gap, /*condense*/false);
            ++closed_gaps;
        }
        INFO("Closed " << closed_gaps << " gaps");
    }
private:
    DECL_LOGGER("MultiGapJoiner");
};

class HybridGapCloser {
public:
    typedef std::function<string (const vector<string>&)> ConsensusF;
private:
    typedef RtSeq Kmer;
    typedef typename GapStorage::gap_info_it gap_info_it;

    DECL_LOGGER("HybridGapCloser");

    Graph& g_;
    const GapStorage& storage_;
    const size_t min_weight_;
    ConsensusF consensus_;
    const size_t long_seq_limit_;
    const size_t max_consensus_reads_;

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
        DEBUG(gap_variants.size() << " gap closing variants, lengths: " << PrintLengths(gap_variants));
        DEBUG("var size original " << gap_variants.size());
        vector<string> new_gap_variants(gap_variants.begin(), gap_variants.end());
        new_gap_variants.resize(std::min(max_consensus_reads_, gap_variants.size()));
        auto s = consensus_(new_gap_variants);
        DEBUG("consenus for " << g_.int_id(start)
                              << " and " << g_.int_id(end)
                              << " found: '" << s << "'");
        return GapDescription(start, end,
                              Sequence(s),
                              edge_gap_start_position, edge_gap_end_position);
    }

    //all gaps guaranteed to correspond to a single edge pair
    GapInfos PadGaps(gap_info_it start, gap_info_it end) const {
        size_t start_min = std::numeric_limits<size_t>::max();
        size_t end_max = 0;
        size_t long_seqs = 0;
        size_t short_seqs = 0;
        for (auto it = start; it != end; ++it) {
            const auto& gap = *it;
            if (gap.gap_seq.size() > long_seq_limit_)
                long_seqs++;
            else
                short_seqs++;

            start_min = std::min(start_min, gap.edge_gap_start_position);
            end_max = std::max(end_max, gap.edge_gap_end_position);
        }

        const bool exclude_long_seqs = (short_seqs >= min_weight_ && short_seqs > long_seqs);

        GapInfos answer;
        for (auto it = start; it != end; ++it) {
            const auto& gap = *it;

            if (exclude_long_seqs && gap.gap_seq.size() > long_seq_limit_)
                continue;

            string s = g_.EdgeNucls(gap.start).Subseq(start_min + g_.k(), gap.edge_gap_start_position + g_.k()).str();
            s += gap.gap_seq.str();
            s += g_.EdgeNucls(gap.end).Subseq(gap.edge_gap_end_position, end_max).str();
            answer.push_back(GapDescription(gap.start, gap.end, Sequence(s), start_min, end_max));
        }
        return answer;
    }

    GapDescription ConstructConsensus(gap_info_it start_it, gap_info_it end_it) const {
        DEBUG("Considering extension " << g_.str(start_it->end));
        size_t cur_len = end_it - start_it;

        //low weight connections filtered earlier
        VERIFY(cur_len >= min_weight_);

        auto padded_gaps = PadGaps(start_it, end_it);
        //all start and end positions are equal here
        if (padded_gaps.size() < min_weight_) {
            DEBUG("Connection weight too low after padding");
            return INVALID_GAP;
        }

        vector<string> gap_variants;
        std::transform(padded_gaps.begin(), padded_gaps.end(), std::back_inserter(gap_variants), 
                       [](const GapDescription& gap) {
            return gap.gap_seq.str();
        });

        //for (auto it = start_it; it != end_it; ++it) {
        //    VERIFY(it->start == start_it->start);
        //    VERIFY(it->end == start_it->end);
        //    VERIFY(it->edge_gap_start_position == start_it->edge_gap_start_position);
        //    VERIFY(it->edge_gap_end_position == start_it->edge_gap_end_position);
        //}
        auto padded_gap = padded_gaps.front();

        return ConstructConsensus(padded_gap.start, padded_gap.end,
                                  padded_gap.edge_gap_start_position,
                                  padded_gap.edge_gap_end_position,
                                  gap_variants);
    }

    GapDescription ConstructConsensus(EdgeId e) const {
        DEBUG("Constructing consensus for edge " << g_.str(e));
        vector<GapDescription> closures;
        for (const auto& edge_pair_gaps : storage_.EdgePairGaps(get(storage_.inner_index(), e))) {
            auto consensus = ConstructConsensus(edge_pair_gaps.first, edge_pair_gaps.second);
            if (consensus != INVALID_GAP) {
                closures.push_back(consensus);
            }
        }

        if (closures.size() == 1) {
            DEBUG("Found unique extension " << closures.front().str(g_));
            return closures.front();
        }

        if (closures.size() > 1) {
            DEBUG("Non-unique extension");
        }
        return INVALID_GAP;
    }

    vector<GapDescription> ConstructConsensus() const {
        vector<vector<GapDescription>> closures_by_thread(omp_get_max_threads());

        # pragma omp parallel for
        for (size_t i = 0; i < storage_.size(); i++) {
            EdgeId e = storage_[i];
            size_t thread_num = omp_get_thread_num();
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

public:
    HybridGapCloser(Graph& g, const GapStorage& storage,
                    size_t min_weight, ConsensusF consensus,
                    size_t long_seq_limit,
                    size_t max_consensus_reads = 20)
            : g_(g), storage_(storage),
              min_weight_(min_weight),
              consensus_(consensus),
              long_seq_limit_(long_seq_limit),
              max_consensus_reads_(max_consensus_reads) {
    }

    map<EdgeId, EdgeId> operator()() {
        EdgeFateTracker fate_tracker(g_);
        MultiGapJoiner gap_joiner(g_);

        gap_joiner(ConstructConsensus());

        CompressAllVertices(g_, true, /*chunk_cnt*/100);
        return fate_tracker.Old2NewMapping();
    };

};

}
}
