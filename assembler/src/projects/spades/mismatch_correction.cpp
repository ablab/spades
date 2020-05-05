//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "mismatch_correction.hpp"

#include "modules/alignment/sequence_mapper.hpp"
#include "modules/simplification/compressor.hpp"

#include "io/reads/read_stream_vector.hpp"
#include "io/dataset_support/read_converter.hpp"

#include "pipeline/graph_pack.hpp"
#include "utils/logger/logger.hpp"

#include <parallel_hashmap/phmap.h>

namespace debruijn_graph {

namespace mismatches {
struct NuclCount {
    std::array<size_t, 4> counts_;

    NuclCount()
            : counts_{} {}

    size_t &operator[](size_t nucl) {
        return counts_[nucl];
    }

    size_t operator[](size_t nucl) const {
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
        for (const auto &entry : other.info_) {
            info_[entry.first] += entry.second;
        }
    }

    void IncIfContains(size_t position, size_t nucl) {
        auto it = info_.find(position);
        if (it != info_.end()) {
            it->second[nucl] += 1;
        }
    }

    void AddPosition(size_t position) {
        info_[position]; //in case map did not contain this key creates entry in the map with default value
    }

public:
    phmap::flat_hash_map<size_t, NuclCount> info_;
};

class MismatchStatistics {
private:
    typedef Graph::EdgeId EdgeId;
    typedef phmap::node_hash_map<EdgeId, MismatchEdgeInfo> InnerMismatchStatistics;
    typedef typename InnerMismatchStatistics::const_iterator const_iterator;
    InnerMismatchStatistics statistics_;

    void CollectPotentialMismatches(const GraphPack &gp) {
        const auto &kmer_mapper = gp.get<KmerMapper<Graph>>();
        const auto &index = gp.get<EdgeIndex<Graph>>();

        for (const auto &mentry : kmer_mapper) {
            // Kmer mapper iterator dereferences to pair (KMer, KMer), not to the reference!
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
                    if (from[i] != to[i] && index.contains(to)) {
                        const auto &position = index.get(to);
                        //FIXME add only canonical edges?
                        statistics_[position.first].AddPosition(position.second + i);
                    }
                }
            }
        }
    }

    void Merge(const MismatchStatistics &other) {
        for (const auto &e_info : other.statistics_) {
            statistics_[e_info.first] += e_info.second;
        }
    }

public:
    MismatchStatistics(const GraphPack &gp) {
        CollectPotentialMismatches(gp);
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

    template<class read_type>
    void Count(io::ReadStream<read_type> &stream, const GraphPack &gp) {
        stream.reset();
        DEBUG("count started");
        auto sm = MapperInstance(gp);
        DEBUG("seq mapper created");
        read_type read;
        const auto &graph = gp.get<Graph>();
        while (!stream.eof()) {
            stream >> read;
            const Sequence &s_read = read.sequence();
            omnigraph::MappingPath<EdgeId> path = sm->MapSequence(s_read,
                                                                  true /* only_simple */);
            TRACE("read mapped");
            VERIFY(path.size() <= 1);
            if (path.size() != 1)
                continue;

            EdgeId e = path[0].first;
            MappingRange mr = path[0].second;

            if (mr.initial_range.size() == mr.mapped_range.size()) {
                const Sequence &s_edge = graph.EdgeNucls(e);
                size_t len = mr.initial_range.size() + graph.k();
                size_t cnt = 0;
                for (size_t i = 0; i < len; i++) {
                    if (s_read[mr.initial_range.start_pos + i] != s_edge[mr.mapped_range.start_pos + i]) {
                        cnt++;
                    }
                }
                if (cnt <= graph.k() / 3) {
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
    void ParallelCount(SingleStreamList &streams, const GraphPack &gp) {
        std::vector<MismatchStatistics> statistics(streams.size(), *this);

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

    GraphPack &gp_;
    Graph &graph_;
    const size_t k_;
    const double relative_threshold_;

    EdgeId CorrectNucl(EdgeId edge, size_t position, char nucl) {
        VERIFY(position >= k_);
        if (position + 1 < graph_.length(edge)) {
            auto tmp = graph_.SplitEdge(edge, position + 1);
            edge = tmp.first;
        }
        EdgeId mismatch = edge;
        if (position > k_) {
            auto tmp = graph_.SplitEdge(edge, position - k_);
            edge = tmp.first;
            mismatch = tmp.second;
        }
        Sequence s_mm = graph_.EdgeNucls(mismatch);
        Sequence correct = s_mm.Subseq(0, k_) + Sequence(std::string(1, nucl)) +
                           s_mm.Subseq(k_ + 1, k_ * 2 + 1);

        VERIFY(nucl != s_mm[k_]);
        EdgeId correct_edge = graph_.AddEdge(graph_.EdgeStart(mismatch), graph_.EdgeEnd(mismatch), correct);
        EdgeId glued = graph_.GlueEdges(mismatch, correct_edge);
        return position > k_ ? edge : glued;
    }

    EdgeId CorrectNucls(EdgeId edge, const std::vector<std::pair<size_t, char>> &mismatches) {
        // Nothing to correct, bail out.
        // Note that this might be a correctness thing as well, as we're calling Compress
        // down below.
        if (mismatches.empty())
            return edge;

        for (auto it = mismatches.rbegin(); it != mismatches.rend(); ++it) {
            edge = CorrectNucl(edge, it->first, it->second);
        }
        EdgeId tmp = omnigraph::Compressor<Graph>(graph_).CompressVertexEdgeId(graph_.EdgeEnd(edge));
        if (tmp == EdgeId())
            return edge;
        else
            return tmp;
    }

    std::vector<std::pair<size_t, char>> FindMismatches(EdgeId edge, const MismatchEdgeInfo &statistics) const {
        std::vector<std::pair<size_t, char>> to_correct;
        const Sequence &s_edge = graph_.EdgeNucls(edge);
        for (size_t i = k_; i < graph_.length(edge); i++) {
            size_t cur_best = 0;
            NuclCount nc = statistics[i];
            for (size_t j = 1; j < 4; j++) {
                if (nc[j] > nc[cur_best]) {
                    cur_best = j;
                }
            }
            char nucl_code = s_edge[i];
            if ((double) nc[cur_best] > relative_threshold_ * (double) nc[nucl_code] + 1.) {
                to_correct.emplace_back(i, cur_best);
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

    size_t CorrectAllEdges(const MismatchStatistics &statistics) {
        size_t res = 0;
        btree::btree_set<EdgeId> conjugate_fix;
        //FIXME after checking saves replace with
        //for (auto it = g_.ConstEdgeBegin(/*canonical only*/true); !it.IsEnd(); ++it) {

#if 0        
        for (EdgeId e : g_.edges()) {
            if (!conjugate_fix.count(g_.conjugate(e)))
                conjugate_fix.insert(e);
        }
#else
        for (auto it = graph_.ConstEdgeBegin(); !it.IsEnd(); ++it) {
            if (conjugate_fix.find(graph_.conjugate(*it)) == conjugate_fix.end()) {
                conjugate_fix.insert(*it);
            }
        }
#endif
        for (EdgeId e : conjugate_fix) {
            DEBUG("processing edge" << graph_.int_id(e));

            auto stat_it = statistics.find(e);
            if (stat_it == statistics.end())
                continue;

            if (!graph_.RelatedVertices(graph_.EdgeStart(e), graph_.EdgeEnd(e))) {
                res += CorrectEdge(e, stat_it->second);
            }
        }
        INFO("All edges processed");
        return res;
    }

    template<class SingleReadStream>
    size_t StopMismatchIteration(SingleReadStream &stream) {
        MismatchStatistics statistics(gp_);
        statistics.Count(stream, gp_);
        return CorrectAllEdges(statistics);
    }

    template<class SingleReadStreamList>
    size_t ParallelStopMismatchIteration(SingleReadStreamList &streams) {
        MismatchStatistics statistics(gp_);
        statistics.ParallelCount(streams, gp_);
        return CorrectAllEdges(statistics);
    }

public:
    MismatchShallNotPass(GraphPack &gp, double relative_threshold = 1.5) :
            gp_(gp),
            graph_(gp.get_mutable<Graph>()),
            k_(gp.k()),
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

} // namespace mismatches

void MismatchCorrection::run(GraphPack &gp, const char*) {
    gp.EnsureBasicMapping();

    auto& dataset = cfg::get_writable().ds;
    std::vector<size_t> libs;
    for (size_t i = 0; i < dataset.reads.lib_count(); ++i) {
        if (dataset.reads[i].is_mismatch_correctable())
            libs.push_back(i);
    }
    auto streams = io::single_binary_readers_for_libs(dataset.reads, libs);
    size_t corrected = mismatches::MismatchShallNotPass(gp, 2).
                       ParallelStopAllMismatches(streams, 1);
    INFO("Corrected " << corrected << " nucleotides");
}

}
