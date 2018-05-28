//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/graph_pack.hpp"
#include "modules/simplification/compressor.hpp"
#include "assembly_graph/handlers/id_track_handler.hpp"
#include "utils/logger/logger.hpp"

#include "io/reads/read_stream_vector.hpp"
#include "modules/alignment/sequence_mapper.hpp"

#include "pipeline/config_struct.hpp"

namespace debruijn_graph {

namespace mismatches {
struct NuclCount {
    std::array<size_t, 4> counts_;

    NuclCount() : counts_{} {
    }

    size_t &operator[](size_t nucl) {
        return counts_[nucl];
    }

    NuclCount &operator+=(const NuclCount &other) {
        counts_[0] += other.counts_[0];
        counts_[1] += other.counts_[1];
        counts_[2] += other.counts_[2];
        counts_[3] += other.counts_[3];
        return *this;
    }
};

struct MismatchEdgeInfo {
    NuclCount operator[](size_t i) const {
        auto it = info_.find(i);
        if (it == info_.end())
            return NuclCount();
        else
            return it->second;
    }

    void operator+=(const MismatchEdgeInfo &other) {
        for (auto it = other.info_.begin(); it != other.info_.end(); ++it) {
            info_[it->first] += it->second;
        }
    }

    void IncIfContains(size_t position, size_t nucl) {
        auto it = info_.find(position);
        if (it != info_.end()) {
            it->second[nucl]++;
        }
    }

    void AddPosition(size_t position) {
        info_[position]; //in case map did not contain this key creates entry in the map with default value
    }

public:
    map<size_t, NuclCount> info_;
};

template<typename EdgeId>
class MismatchStatistics {
private:
    typedef typename map<EdgeId, MismatchEdgeInfo>::const_iterator const_iterator;
    map<EdgeId, MismatchEdgeInfo> statistics_;

    template<class graph_pack>
    void CollectPotensialMismatches(const graph_pack &gp) {
        const auto &kmer_mapper = gp.kmer_mapper;
        for (auto it = kmer_mapper.begin(); it != kmer_mapper.end(); ++it) {
            // Kmer mapper iterator dereferences to pair (KMer, KMer), not to the reference!
            const auto mentry = *it;
            const RtSeq &from = mentry.first;
            const RtSeq &to = mentry.second;
            size_t cnt = 0;
            std::array<size_t, 4> cnt_arr{};

            for (size_t i = 0; i < from.size(); i++) {
                if (from[i] != to[i]) {
                    cnt++;
                    cnt_arr[(i * 4) / from.size()]++;
                }
            }

            //last two conditions - to avoid excessive indels.
            //if two/third of nucleotides in first/last quarter are mismatches, then it means erroneous mapping
            if (cnt >= 1 && cnt <= from.size() / 3 && cnt_arr[0] <= from.size() / 6 &&
                cnt_arr[3] <= from.size() / 6) {
                for (size_t i = 0; i < from.size(); i++) {
                    if (from[i] != to[i] && gp.index.contains(to)) {
                        pair<EdgeId, size_t> position = gp.index.get(to);
                        //FIXME add only canonical edges?
                        statistics_[position.first].AddPosition(position.second + i);
                    }
                }
            }
        }
    }

    void Merge(const MismatchStatistics<EdgeId> &other) {
        for (const auto &e_info : other.statistics_) {
            statistics_[e_info.first] += e_info.second;
        }
    }

public:
    MismatchStatistics(const conj_graph_pack &gp) {
        CollectPotensialMismatches(gp);
    }

    const_iterator begin() const {
        return statistics_.begin();
    }

    const_iterator end() const {
        return statistics_.end();
    }

    const_iterator find(const EdgeId &edge) const {
        return statistics_.find(edge);
    }

    template<class graph_pack, class read_type>
    void Count(io::ReadStream<read_type> &stream, const graph_pack &gp) {
        stream.reset();
        DEBUG("count started");
        auto sm = MapperInstance(gp);
        DEBUG("seq mapper created");
        read_type read;
        while (!stream.eof()) {
            stream >> read;
            const Sequence &s_read = read.sequence();
            omnigraph::MappingPath<EdgeId> path = sm->MapSequence(s_read);
            TRACE("read mapped");
            if (path.size() != 1)
                continue;

            EdgeId e = path[0].first;
            MappingRange mr = path[0].second;

            if (mr.initial_range.size() == mr.mapped_range.size()) {
                const Sequence &s_edge = gp.g.EdgeNucls(e);
                size_t len = mr.initial_range.size() + gp.g.k();
                size_t cnt = 0;
                for (size_t i = 0; i < len; i++) {
                    if (s_read[mr.initial_range.start_pos + i] != s_edge[mr.mapped_range.start_pos + i]) {
                        cnt++;
                    }
                }
                if (cnt <= gp.g.k() / 3) {
                    TRACE("statistics changing");
                    auto it = statistics_.find(e);
                    if (it == statistics_.end()) {
                        //                            if (gp.g.length(path[0].first) < 4000)
                        //                                WARN ("id "<< gp.g.length(path[0].first)<<"  " << len);
                        continue;
                    }
                    for (size_t i = 0; i < len; i++) {
                        char nucl_code = s_read[mr.initial_range.start_pos + i];
                        it->second.IncIfContains(mr.mapped_range.start_pos + i, nucl_code);
                    }
                }
            }
        }
    }

    template<class SingleStreamList>
    void ParallelCount(SingleStreamList &streams, const conj_graph_pack &gp) {
        std::vector<MismatchStatistics<EdgeId>> statistics(streams.size(), *this);

        #pragma omp parallel for
        for (size_t i = 0; i < streams.size(); ++i) {
            DEBUG("statistics created thread " << i);
            statistics[i].Count(streams[i], gp);
            DEBUG("count finished thread " << i);
        }

        INFO("Finished collecting potential mismatches positions");
        for (size_t i = 0; i < statistics.size(); i++) {
            Merge(statistics[i]);
        }
    }
};

class MismatchShallNotPass {
private:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    conj_graph_pack &gp_;
    Graph &g_;
    const size_t k_;
    const double relative_threshold_;

    EdgeId CorrectNucl(EdgeId edge, size_t position, char nucl) {
        VERIFY(position >= k_);
        if (position + 1 < g_.length(edge)) {
            edge = g_.SplitEdge(edge, position + 1).first;
        }
        EdgeId mismatch = edge;
        if (position > k_) {
            auto tmp = g_.SplitEdge(edge, position - k_);
            edge = tmp.first;
            mismatch = tmp.second;
        }
        Sequence s_mm = g_.EdgeNucls(mismatch);
        Sequence correct = s_mm.Subseq(0, k_) + Sequence(string(1, nucl)) +
                           s_mm.Subseq(k_ + 1, k_ * 2 + 1);

        VERIFY(nucl != s_mm[k_]);
        EdgeId correct_edge = g_.AddEdge(g_.EdgeStart(mismatch), g_.EdgeEnd(mismatch), correct);
        EdgeId glued = g_.GlueEdges(mismatch, correct_edge);
        return position > k_ ? edge : glued;
    }

    EdgeId CorrectNucls(EdgeId edge, const std::vector<pair<size_t, char>> &mismatches) {
        for (auto it = mismatches.rbegin(); it != mismatches.rend(); ++it) {
            edge = CorrectNucl(edge, it->first, it->second);
        }
        EdgeId tmp = Compressor<Graph>(g_).CompressVertexEdgeId(g_.EdgeEnd(edge));
        if (tmp == EdgeId())
            return edge;
        else
            return tmp;
    }

    vector<pair<size_t, char>> FindMismatches(EdgeId edge, const MismatchEdgeInfo &statistics) const {
        vector<pair<size_t, char>> to_correct;
        const Sequence &s_edge = g_.EdgeNucls(edge);
        for (size_t i = k_; i < g_.length(edge); i++) {
            size_t cur_best = 0;
            NuclCount nc = statistics[i];
            for (size_t j = 1; j < 4; j++) {
                if (nc[j] > nc[cur_best]) {
                    cur_best = j;
                }
            }
            char nucl_code = s_edge[i];
            if ((double) nc[cur_best] > relative_threshold_ * (double) nc[nucl_code] + 1.) {
                to_correct.push_back(make_pair(i, cur_best));
                i += k_;
            }

        }
        return to_correct;
    }

    size_t CorrectEdge(EdgeId edge, const MismatchEdgeInfo &statistics) {
        auto to_correct = FindMismatches(edge, statistics);
        EdgeId new_edge = CorrectNucls(edge, to_correct);
        if (new_edge == EdgeId())
            new_edge = edge;

        return to_correct.size();
    }

    size_t CorrectAllEdges(const MismatchStatistics<EdgeId> &statistics) {
        size_t res = 0;
        set<EdgeId> conjugate_fix;
        //FIXME after checking saves replace with
        //for (auto it = g_.ConstEdgeBegin(/*canonical only*/true); !it.IsEnd(); ++it) {

        for (auto it = gp_.g.ConstEdgeBegin(); !it.IsEnd(); ++it) {
            if (conjugate_fix.find(gp_.g.conjugate(*it)) == conjugate_fix.end()) {
                conjugate_fix.insert(*it);
            }
        }
        for (auto it = conjugate_fix.begin(); it != conjugate_fix.end(); ++it) {
            EdgeId e = *it;
            DEBUG("processing edge" << g_.int_id(e));

            auto stat_it = statistics.find(e);
            if (stat_it != statistics.end()) {
                if (!g_.RelatedVertices(g_.EdgeStart(e), g_.EdgeEnd(e))) {
                    res += CorrectEdge(e, stat_it->second);
                }
            }
        }
        INFO("All edges processed");
        return res;
    }

    template<class SingleReadStream>
    size_t StopMismatchIteration(SingleReadStream &stream) {
        MismatchStatistics<EdgeId> statistics(gp_);
        statistics.Count(stream, gp_);
        return CorrectAllEdges(statistics);
    }

    template<class SingleReadStreamList>
    size_t ParallelStopMismatchIteration(SingleReadStreamList &streams) {
        MismatchStatistics<EdgeId> statistics(gp_);
        statistics.ParallelCount(streams, gp_);
        return CorrectAllEdges(statistics);
    }

public:
    MismatchShallNotPass(conj_graph_pack &gp, double relative_threshold = 1.5) :
            gp_(gp),
            g_(gp.g),
            k_(gp.g.k()),
            relative_threshold_(relative_threshold) {
        VERIFY(relative_threshold >= 1);
    }


    template<class SingleReadStream>
    size_t StopAllMismatches(SingleReadStream &stream, size_t max_iterations = 1) {
        size_t res = 0;
        while (max_iterations > 0) {
            size_t last = StopMismatchIteration(stream);
            res += last;
            if (last == 0)
                break;
            max_iterations--;
        }
        return res;
    }

    template<class SingleReadStreamList>
    size_t ParallelStopAllMismatches(SingleReadStreamList &streams, size_t max_iterations = 1) {
        size_t res = 0;
        while (max_iterations > 0) {
            size_t last = ParallelStopMismatchIteration(streams);
            res += last;
            if (last == 0)
                break;
            max_iterations--;
        }
        return res;
    }
};

}
}
