//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "modules/simplification/compressor.hpp"
#include "assembly_graph/handlers/id_track_handler.hpp"
#include "utils/logger/logger.hpp"

#include "io/reads/read_stream_vector.hpp"
#include "modules/alignment/sequence_mapper.hpp"

#include "pipeline/config_struct.hpp"

namespace debruijn_graph {

namespace mismatches {
struct NuclCount {
    size_t counts_[4];

    NuclCount() {
        memset(counts_, 0, sizeof(counts_));
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
        auto &kmer_mapper = gp.kmer_mapper;
        for (auto it = kmer_mapper.begin(); it != kmer_mapper.end(); ++it) {
            // Kmer mapper iterator dereferences to pair (KMer, KMer), not to the reference!
            const auto mentry = *it;
            const RtSeq &from = mentry.first;
            const RtSeq &to = mentry.second;
            size_t cnt = 0;
            size_t cnt_arr[4];
            for (size_t i = 0; i < 4; i++)
                cnt_arr[i] = 0;
            for (size_t i = 0; i < from.size(); i++) {
                if (from[i] != to[i]) {
                    cnt++;
                    cnt_arr[(i * 4) / from.size()]++;
                }
            }
            //last two contitions - to avoid excessive indels.
            //if two/third of nucleotides in first/last quoter are mismatches, then it means erroneous mapping

            if (cnt >= 1 && cnt <= from.size() / 3 && cnt_arr[0] <= from.size() / 6 &&
                cnt_arr[3] <= from.size() / 6) {
                for (size_t i = 0; i < from.size(); i++) {
                    if (from[i] != to[i] && gp.index.contains(to)) {
                        pair<EdgeId, size_t> position = gp.index.get(to);
                        statistics_[position.first].AddPosition(position.second + i);
                    }
                }
            }
        }
    }

    void operator+=(const MismatchStatistics<EdgeId> &other) {
        for (auto it = other.statistics_.begin(); it != other.statistics_.end(); ++it) {
            statistics_[it->first] += it->second;
        }
    }

public:
    template<class graph_pack>
    MismatchStatistics(const graph_pack &gp) {
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
        while (!stream.eof()) {
            read_type read;
            stream >> read;
            const Sequence &s_read = read.sequence();
            omnigraph::MappingPath<EdgeId> path = sm->MapSequence(s_read);
            TRACE("read mapped");
            if (path.size() == 1 && path[0].second.initial_range.size() == path[0].second.mapped_range.size()) {
                Range initial_range = path[0].second.initial_range;
                Range mapped_range = path[0].second.mapped_range;
                const Sequence &s_edge = gp.g.EdgeNucls(path[0].first);
                size_t len = initial_range.size() + gp.g.k();
                size_t cnt = 0;
                for (size_t i = 0; i < len; i++) {
                    if (s_read[initial_range.start_pos + i] != s_edge[mapped_range.start_pos + i]) {
                        cnt++;
                    }
                }
                if (cnt <= gp.g.k() / 3) {
                    TRACE("statistics changing");
                    auto it = statistics_.find(path[0].first);
                    if (it == statistics_.end()) {
                        //                            if (gp.g.length(path[0].first) < 4000)
                        //                                WARN ("id "<< gp.g.length(path[0].first)<<"  " << len);
                        continue;
                    }
                    for (size_t i = 0; i < len; i++) {
                        size_t nucl_code = s_read[initial_range.start_pos + i];
                        it->second.IncIfContains(mapped_range.start_pos + i, nucl_code);
                    }
                }
            }
        }
    }

    template<class graph_pack, class read_type>
    void ParallelCount(io::ReadStreamList<read_type> &streams, const graph_pack &gp) {
        size_t nthreads = streams.size();
        std::vector<MismatchStatistics<EdgeId> *> statistics(nthreads);
#pragma omp parallel for num_threads(nthreads) shared(streams, statistics)
        for (size_t i = 0; i < nthreads; ++i) {
            statistics[i] = new MismatchStatistics<EdgeId>(*this);
            DEBUG("statistics created thread " << i);
            statistics[i]->Count(streams[i], gp);
            DEBUG("count finished thread " << i);
        }

        INFO("Finished collecting potential mismatches positions");
        for (size_t i = 0; i < statistics.size(); i++) {
            *this += *statistics[i];
            delete statistics[i];
        }
    }
};
}

template<class graph_pack, class read_type>
class MismatchShallNotPass {
private:
    typedef typename graph_pack::graph_t Graph;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    graph_pack &gp_;
    double relative_threshold_;

    EdgeId CorrectNucl(EdgeId edge, size_t position, char nucl) {
        VERIFY(position >= gp_.g.k());
        if (position + 1 < gp_.g.length(edge)) {
            edge = gp_.g.SplitEdge(edge, position + 1).first;
        }
        EdgeId mismatch = edge;
        if (position > gp_.g.k()) {
            auto tmp = gp_.g.SplitEdge(edge, position - gp_.g.k());
            edge = tmp.first;
            mismatch = tmp.second;
        }
        Sequence s_mm = gp_.g.EdgeNucls(mismatch);
        Sequence correct = s_mm.Subseq(0, gp_.g.k()) + Sequence(string(1, nucl)) +
                           s_mm.Subseq(gp_.g.k() + 1, gp_.g.k() * 2 + 1);

        VERIFY(nucl != s_mm[gp_.g.k()]);
        EdgeId correct_edge = gp_.g.AddEdge(gp_.g.EdgeStart(mismatch), gp_.g.EdgeEnd(mismatch), correct);
        EdgeId glued = gp_.g.GlueEdges(mismatch, correct_edge);
        return position > gp_.g.k() ? edge : glued;
    }

    EdgeId CorrectNucls(EdgeId edge, const std::vector<pair<size_t, char>> &mismatches) {
        for (auto it = mismatches.rbegin(); it != mismatches.rend(); ++it) {
            edge = CorrectNucl(edge, it->first, it->second);
        }
        EdgeId tmp = Compressor<Graph>(gp_.g).CompressVertexEdgeId(gp_.g.EdgeEnd(edge));
        if (tmp == EdgeId(0))
            return edge;
        else
            return tmp;
    }

    vector<pair<size_t, char>> FindMismatches(EdgeId edge, const mismatches::MismatchEdgeInfo &statistics) {
        vector<pair<size_t, char>> to_correct;
        const Sequence &s_edge = gp_.g.EdgeNucls(edge);
        for (size_t i = gp_.g.k(); i < gp_.g.length(edge); i++) {
            size_t cur_best = 0;
            mismatches::NuclCount nc = statistics[i];
            for (size_t j = 1; j < 4; j++) {
                if (nc[j] > nc[cur_best]) {
                    cur_best = j;
                }
            }
            size_t nucl_code = s_edge[i];
            if ((double) nc[cur_best] > relative_threshold_ * (double) nc[nucl_code] + 1.) {
                to_correct.push_back(make_pair(i, cur_best));
                i += gp_.g.k();
            }

        }
        return to_correct;
    }

    size_t CorrectEdge(EdgeId edge, const mismatches::MismatchEdgeInfo &statistics) {
        vector<pair<size_t, char>> to_correct = FindMismatches(edge, statistics);
        EdgeId new_edge = CorrectNucls(edge, to_correct);
        if (new_edge == EdgeId(0))
            new_edge = edge;

        return to_correct.size();
    }

    size_t CorrectAllEdges(const mismatches::MismatchStatistics<typename Graph::EdgeId> &statistics) {
        size_t res = 0;
        set<EdgeId> conjugate_fix;
        for (auto it = gp_.g.ConstEdgeBegin(); !it.IsEnd(); ++it) {
            if (conjugate_fix.find(gp_.g.conjugate(*it)) == conjugate_fix.end()) {
                conjugate_fix.insert(*it);
            }
        }
        for (auto it = conjugate_fix.begin(); it != conjugate_fix.end(); ++it) {
            DEBUG("processing edge" << gp_.g.int_id(*it));

            if (statistics.find(*it) != statistics.end()) {
                if (!gp_.g.RelatedVertices(gp_.g.EdgeStart(*it), gp_.g.EdgeEnd(*it)))
                    res += CorrectEdge(*it, statistics.find(*it)->second);
            }
        }
        INFO("All edges processed");
        return res;
    }

    size_t StopMismatchIteration(io::ReadStream<read_type> &stream) {
        mismatches::MismatchStatistics<typename Graph::EdgeId> statistics(gp_);
        statistics.Count(stream, gp_);
        return CorrectAllEdges(statistics);
    }

    size_t ParallelStopMismatchIteration(io::ReadStreamList<read_type> &streams) {
        mismatches::MismatchStatistics<typename Graph::EdgeId> statistics(gp_);
        statistics.ParallelCount(streams, gp_);
        return CorrectAllEdges(statistics);
    }

public:
    MismatchShallNotPass(graph_pack &gp, double relative_threshold = 1.5) :
            gp_(gp),
            relative_threshold_(relative_threshold) {
        VERIFY(relative_threshold >= 1);
    }


    size_t StopAllMismatches(io::ReadStream<read_type> &stream, size_t max_iterations = 1) {
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

    size_t ParallelStopAllMismatches(io::ReadStreamList<read_type> &streams, size_t max_iterations = 1) {
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
