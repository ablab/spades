//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "mismatch_correction.hpp"

#include "modules/alignment/sequence_mapper.hpp"
#include "modules/alignment/sequence_mapper_notifier.hpp"
#include "modules/simplification/compressor.hpp"

#include "io/reads/read_stream_vector.hpp"
#include "io/dataset_support/read_converter.hpp"

#include "pipeline/graph_pack.hpp"
#include "utils/logger/logger.hpp"

#include "adt/flat_set.hpp"
#include <parallel_hashmap/phmap.h>

template <typename Iter>
std::vector<Iter> split_iterator(size_t chunks, Iter b, Iter e, size_t n) {
    std::vector<Iter> result(chunks + 1, e);
    size_t leading_chunks_size = n / chunks;
    if (n % chunks) {
        leading_chunks_size += 1;
    }
    DEBUG("Leading chunk size " << leading_chunks_size << " n " << n);
    for (size_t i = 0; b != e; ++b, ++i) {
        if (i % leading_chunks_size == 0) {
            DEBUG("i = " << i << ", setting to " << i / leading_chunks_size);
            result[i / leading_chunks_size] = b;
            VERIFY(i / leading_chunks_size < chunks);
        }
    }
    result[chunks] = e;

    return result;
}


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
        auto it = info_.find(uint32_t(i));
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

    void increment(size_t position, size_t nucl) {
        info_[uint32_t(position)][nucl] += 1;
    }

    void ClearValues() {
        for (auto &kv : info_) {
            kv.second = NuclCount();
        }
    }

public:
    phmap::flat_hash_map<uint32_t, NuclCount> info_;
};

class MismatchStatistics : public SequenceMapperListener {
private:
    typedef Graph::EdgeId EdgeId;
    typedef phmap::node_hash_map<EdgeId, MismatchEdgeInfo> InnerMismatchStatistics;
    typedef typename InnerMismatchStatistics::const_iterator const_iterator;
    InnerMismatchStatistics statistics_;
    std::vector<InnerMismatchStatistics> statistics_buffers_;

    typedef phmap::node_hash_map<EdgeId, adt::flat_set<uint32_t>> MismatchCandidates;
    MismatchCandidates candidates_;

    const Graph &g_;

    template <typename Iter>
    void CollectPotentialMismatches(const GraphPack &gp, Iter b, Iter e, MismatchCandidates &candidates) {
        const auto &index = gp.get<EdgeIndex<Graph>>();

        for (const auto &mentry : adt::make_range(b, e)) {
            // Kmer mapper iterator dereferences to pair (KMer, KMer), not to the reference!
            const RtSeq &from = mentry.first;
            const RtSeq &to = mentry.second;

            // No need to do anything if the target is not in the graph.
            // This certainly expects normalized index.
            if (!index.contains(to))
                continue;

            size_t cnt = 0;
            std::array<size_t, 4> cnt_arr{};

            for (size_t i = 0; i < from.size(); i++) {
                if (from[i] == to[i])
                    continue;

                cnt += 1;
                cnt_arr[(i * 4) / from.size()] += 1;
            }

            // No mismatches, no cookies
            if (cnt == 0)
                continue;

            // If there are too many mismatches, then it means erroneous mapping
            if (cnt > from.size() / 3)
                continue;

            // These conditions are to avoid excessive indels: if two/third of
            // nucleotides in first/last quarter are mismatches, then it means
            // erroneous mapping
            if (cnt_arr[0] > from.size() / 6 ||  cnt_arr[3] > from.size() / 6)
                continue;

            const auto &position = index.get(to);
            for (size_t i = 0; i < from.size(); i++) {
                if (from[i] == to[i])
                    continue;

                if (position.second > std::numeric_limits<uint32_t>::max())
                    continue;

                //FIXME add only canonical edges?
                candidates[position.first].insert(uint32_t(position.second + i));
            }
        }
    }

    void CollectPotentialMismatches(const GraphPack &gp) {
        size_t nthreads = omp_get_max_threads();
        const auto &kmer_mapper = gp.get<KmerMapper<Graph>>();
        auto iters = split_iterator(nthreads, kmer_mapper.begin(), kmer_mapper.end(), kmer_mapper.size());
        VERIFY(iters.front() == kmer_mapper.begin());
        VERIFY(iters.back() == kmer_mapper.end());

        std::vector<MismatchCandidates> potential_mismatches(nthreads);
#       pragma omp parallel for
        for (size_t i = 0; i < nthreads; ++i) {
            CollectPotentialMismatches(gp, iters[i], iters[i + 1], potential_mismatches[i]);
        }

        for (auto &entry : potential_mismatches) {
            for (const auto &candidate : entry) {
                candidates_[candidate.first].insert(candidate.second.begin(),
                                                    candidate.second.end());
            }
            entry.clear();
        }

        {
            size_t edges = candidates_.size();
            size_t positions = 0;
            for (const auto &candidate : candidates_) {
                positions += candidate.second.size();
            }

            INFO("Total " << edges << " edges (out of " << gp.get<Graph>().e_size() <<  ") with " << positions << " potential mismatch positions ("
                 << double(positions) / double(edges) << " positions per edge)");
        }
    }

    template <typename Read>
    void ProcessSingleReadImpl(size_t thread_index, const Read& read, const MappingPath<EdgeId> &path) {
        // VERIFY(path.size() <= 1);
        if (path.size() != 1)  // TODO Use only_simple feature
            return;

        EdgeId e = path[0].first;
        MappingRange mr = path[0].second;
        const Sequence &s_read = read.sequence();
        auto &buffer = statistics_buffers_[thread_index];

        if (mr.initial_range.size() != mr.mapped_range.size())
            return;

        auto it = candidates_.find(e);
        if (it == candidates_.end())
            return;

        const Sequence &s_edge = g_.EdgeNucls(e);
        size_t len = mr.initial_range.size() + g_.k();
        size_t cnt = 0;
        for (size_t i = 0; i < len; i++) {
            cnt += (s_read[mr.initial_range.start_pos + i] !=
                    s_edge[mr.mapped_range.start_pos + i]);
        }

        if (cnt > g_.k() / 3)
            return;

        TRACE("statistics might be changing");
        for (size_t i = 0; i < len; i++) {
            size_t pos = mr.mapped_range.start_pos + i;
            if (pos > std::numeric_limits<uint32_t>::max())
                continue;

            if (!it->second.count(uint32_t(pos)))
                continue;

            char nucl_code = s_read[mr.initial_range.start_pos + i];
            buffer[e].increment(pos, nucl_code);
        }
    }

    void Merge(InnerMismatchStatistics &other_statistics) {
        for (auto &e_info : other_statistics) {
            statistics_[e_info.first] += e_info.second;
            e_info.second.ClearValues();
        }
    }

public:
    MismatchStatistics(const GraphPack &gp):
            g_(gp.get<Graph>()) {
        CollectPotentialMismatches(gp);
    }

    void StartProcessLibrary(size_t threads_count) override {
        statistics_buffers_.clear();
        statistics_buffers_.resize(threads_count, statistics_);
    }

    void StopProcessLibrary() override {
        statistics_buffers_.clear();
    }

    void ProcessSingleRead(size_t thread_index, const io::SingleReadSeq &read, const MappingPath<EdgeId> &path) override {
        ProcessSingleReadImpl(thread_index, read, path);
    }

    void ProcessSingleRead(size_t thread_index, const io::SingleRead &read, const MappingPath<EdgeId> &path) override {
        ProcessSingleReadImpl(thread_index, read, path);
    }

    void MergeBuffer(size_t thread_index) override {
        Merge(statistics_buffers_[thread_index]);
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

    size_t ParallelStopMismatchIteration() {
        INFO("Collect potential mismatches");
        MismatchStatistics statistics(gp_);
        INFO("Potential mismatches collected");

        SequenceMapperNotifier notifier(gp_, cfg::get().ds.reads.lib_count());
        auto& dataset = cfg::get_writable().ds;
        auto mapper = MapperInstance(gp_);

        for (size_t i = 0; i < dataset.reads.lib_count(); ++i) {
            if (!dataset.reads[i].is_mismatch_correctable())
                continue;

            notifier.Subscribe(i, &statistics);
            auto &reads = cfg::get_writable().ds.reads[i];
            auto single_streams = single_binary_readers(reads, /*followed by rc */true, /*binary*/true);
            notifier.ProcessLibrary(single_streams, i, *mapper);
        }

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


    size_t ParallelStopAllMismatches(size_t max_iterations = 1) {
        size_t res = 0;
        while (max_iterations > 0) {
            size_t last = ParallelStopMismatchIteration();
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
    size_t corrected = mismatches::MismatchShallNotPass(gp, 2).
                       ParallelStopAllMismatches(1);
    INFO("Corrected " << corrected << " nucleotides");
}

}
