#pragma once
//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/core/construction_helper.hpp"
#include "kmer_index/extension_index/kmer_extension_index.hpp"
#include "utils/parallel/openmp_wrapper.h"
#include "utils/parallel/parallel_wrapper.hpp"
#include <numeric>

namespace debruijn_graph {

/*
 * Constructs DeBruijnGraph from DeBruijn Graph using "new DeBruijnGraphConstructor(DeBruijn).ConstructGraph(DeBruijnGraph, Index)"
 */
template<class Graph, class Index>
class DeBruijnGraphConstructor {
private:
    typedef typename Graph::EdgeId EdgeId;
    typedef Index DeBruijn;
    typedef typename Graph::VertexId VertexId;
    typedef typename Index::KMer Kmer;
    typedef typename DeBruijn::KeyWithHash KeyWithHash;
    typedef typename DeBruijn::kmer_iterator kmer_iterator;

    Graph &graph_;
    DeBruijn &origin_;
    size_t kmer_size_;

    bool StepRightIfPossible(KeyWithHash &kwh) {
        // VERIFY(origin_.contains(edge));
        if (origin_.RivalEdgeCount(kwh) == 1
                && origin_.NextEdgeCount(kwh) == 1) {
            kwh = origin_.NextEdge(kwh);
            // VERIFY(origin_.contains(next_edge));
            return true;
        }
        return false;
    }

    KeyWithHash &GoRight(KeyWithHash &kwh) {
        KeyWithHash initial = kwh;
        while (StepRightIfPossible(kwh) && kwh != initial) {
            ;
        }
        return kwh;
    }

    KeyWithHash &GoLeft(KeyWithHash &kwh) {
        //These strange things are in order to avoid making copies of kwh
        kwh = !kwh;
        kwh = !GoRight(kwh);
        return kwh;
    }

    Sequence ConstructSeqGoingRight(KeyWithHash &kwh) {
        SequenceBuilder s;
        s.append(kwh.key());
        KeyWithHash initial = kwh;
        while (StepRightIfPossible(kwh) && kwh != initial) {
            s.append(kwh[kmer_size_]);
        }
        return s.BuildSequence();
    }

    Sequence ConstructSequenceWithEdge(const KeyWithHash &kwh) {
        KeyWithHash tmp = kwh;
        return ConstructSeqGoingRight(GoLeft(tmp));
    }

    VertexId FindVertexByOutgoingEdges(Kmer kmer) {
        for (char c = 0; c < 4; ++c) {
            KeyWithHash edge = origin_.ConstructKWH(kmer.pushBack(c));
            if (origin_.contains(edge))
                return graph_.EdgeStart(origin_.get_value(edge).edge_id);
        }
        return VertexId(NULL);
    }

    VertexId FindVertexByIncomingEdges(Kmer kmer) {
        for (char c = 0; c < 4; ++c) {
            KeyWithHash edge = origin_.ConstructKWH(kmer.pushFront(c));
            if (origin_.contains(edge)) {
                return graph_.EdgeEnd(origin_.get_value(edge).edge_id);
            }
        }
        return VertexId(NULL);
    }

    VertexId FindVertex(Kmer kmer) {
        VertexId v = FindVertexByOutgoingEdges(kmer);
        return v == VertexId(NULL) ? FindVertexByIncomingEdges(kmer) : v;
    }

    VertexId FindVertexMaybeMissing(Kmer kmer) {
        VertexId v = FindVertex(kmer);
        return v != VertexId(NULL) ? v : graph_.AddVertex();
    }

    VertexId FindEndMaybeMissing(const ConjugateDeBruijnGraph& graph,
            VertexId start, Kmer start_kmer, Kmer end_kmer) {
        if (start_kmer == end_kmer) {
            return start;
        } else if (start_kmer == !end_kmer) {
            return graph.conjugate(start);
        } else {
            return FindVertexMaybeMissing(end_kmer);
        }
    }

    void ConstructPart(const std::vector<KeyWithHash>& kwh_list,
            std::vector<Sequence>& sequences) {
        for (size_t i = 0; i < sequences.size(); ++i) {
            if (origin_.contains(kwh_list[i])) {
                continue;
            }

            Kmer start_kmer = sequences[i].start < Kmer > (kmer_size_);
            Kmer end_kmer = sequences[i].end < Kmer > (kmer_size_);

            VertexId start = FindVertexMaybeMissing(start_kmer);
            VertexId end = FindEndMaybeMissing(graph_, start, start_kmer,
                    end_kmer);

            graph_.AddEdge(start, end, sequences[i]);
        }
    }

    void AddKmers(kmer_iterator &it, kmer_iterator &end, size_t queueSize,
                  std::vector<KeyWithHash>& kwh_list) {
        for (; kwh_list.size() != queueSize && it != end; ++it) {
            KeyWithHash kwh = origin_.ConstructKWH(Kmer(unsigned(kmer_size_ + 1), (*it).data()));

            if (!origin_.contains(kwh))
                kwh_list.push_back(kwh);
        }
    }

    void CalculateSequences(std::vector<KeyWithHash> &kwh_list,
                            std::vector<Sequence> &sequences) {
        size_t size = kwh_list.size();
        sequences.resize(size);

#       pragma omp parallel for schedule(guided)
        for (size_t i = 0; i < size; ++i) {
            sequences[i] = ConstructSequenceWithEdge(kwh_list[i]);
        }
    }

public:
    DeBruijnGraphConstructor(Graph& graph, DeBruijn &origin) :
            graph_(graph), origin_(origin), kmer_size_(graph_.k()) {
    }

    void ConstructGraph(size_t queueMinSize, size_t queueMaxSize,
                        double queueGrowthRate) {
        kmer_iterator it = origin_.kmer_begin();
        kmer_iterator end = origin_.kmer_end();
        size_t queueSize = queueMinSize;
        std::vector<KeyWithHash> kwh_list;
        std::vector<Sequence> sequences;
        kwh_list.reserve(queueSize);
        sequences.reserve(queueMaxSize);
        while (it != end) {
            AddKmers(it, end, queueSize, kwh_list); // format a queue of kmers that are not in index
            CalculateSequences(kwh_list, sequences); // in parallel
            ConstructPart(kwh_list, sequences);
            kwh_list.clear();
            queueSize = std::min(size_t(double(queueSize) * queueGrowthRate), queueMaxSize);
        }
    }

private:
    DECL_LOGGER("DeBruijnGraphConstructor")
};

class UnbranchingPathExtractor {
private:
    typedef kmers::DeBruijnExtensionIndex<> Index;
    typedef RtSeq Kmer;
    typedef Index::kmer_iterator kmer_iterator;
    typedef Index::DeEdge DeEdge;
    typedef Index::KeyWithHash KeyWithHash;

    Index &origin_;
    size_t kmer_size_;

    bool IsJunction(KeyWithHash kwh) const {
        return IsJunction(origin_.get_value(kwh));
    }

    bool IsJunction(kmers::InOutMask mask) const {
        return !mask.CheckUniqueOutgoing() || !mask.CheckUniqueIncoming();
    }

    void AddStartDeEdgesForVertex(KeyWithHash kh, kmers::InOutMask mask,
                                  std::vector<DeEdge>& start_edges) const {
        for (char next = 0; next < 4; next++) {
            if (!mask.CheckOutgoing(next))
                continue;

            start_edges.emplace_back(kh, origin_.GetOutgoing(kh, next));
            TRACE("Added to queue " << start_edges.back() << " " << mask);
        }
    }

    void AddStartDeEdges(KeyWithHash kh, std::vector<DeEdge>& start_edges) const {
        start_edges.clear();
        auto extensions = origin_.get_value(kh);
        if (!IsJunction(extensions))
            return;

        AddStartDeEdgesForVertex(kh, extensions, start_edges);
        KeyWithHash kh_inv = !kh;
        if (!kh_inv.is_minimal()) {
            AddStartDeEdgesForVertex(kh_inv, origin_.get_value(kh_inv),
                                     start_edges);
        }
    }

    bool StepRightIfPossible(KeyWithHash &kwh) const {
        kmers::InOutMask mask = origin_.get_value(kwh);
        if (mask.CheckUniqueOutgoing() && mask.CheckUniqueIncoming()) {
            kwh = origin_.GetOutgoing(kwh, mask.GetUniqueOutgoing());
            return true;
        }
        return false;
    }

    bool StepRightIfPossible(DeEdge &edge) const {
        kmers::InOutMask mask = origin_.get_value(edge.end);
        if (mask.CheckUniqueOutgoing() && mask.CheckUniqueIncoming()) {
            edge = DeEdge(edge.end,
                          origin_.GetOutgoing(edge.end, mask.GetUniqueOutgoing()));
            return true;
        }
        return false;
    }

    // KeyWithHash has no operator<
    static const KeyWithHash& min_by_key(const KeyWithHash &k1, const KeyWithHash &k2) {
        return k1.key() < k2.key() ? k1 : k2;
    }

    KeyWithHash FindMinimalKMerInLoop(const KeyWithHash &start) const {
        KeyWithHash minimal = min_by_key(start, !start);
        KeyWithHash kh = start;
        StepRightIfPossible(kh);

        for (; start != kh; StepRightIfPossible(kh)) {
            minimal = min_by_key(minimal, kh);
            minimal = min_by_key(minimal, !kh);
        }
        return minimal;
    }

    Sequence ConstructSequenceWithEdge(DeEdge edge, SequenceBuilder &builder) const {
        builder.clear(); // We reuse the buffer to reduce malloc traffic
        builder.append(edge.start.key());
        builder.append(edge.end[kmer_size_ - 1]);
        DeEdge initial = edge;
        while (StepRightIfPossible(edge) && edge != initial) {
            builder.append(edge.end[kmer_size_ - 1]);
        }
        return builder.BuildSequence();
    }

    // Loop consists of 4 parts: 2 selfRC k+1-mers and two sequences of arbitrary length RC to each other; pos is a position of one of selfRC edges
    std::vector<Sequence> SplitLoop(const Sequence &s, size_t pos) const {
        return { s.Subseq(pos, pos + kmer_size_ + 1),
                 s.Subseq(pos + 1, s.size() - kmer_size_) + s.Subseq(0, pos + kmer_size_) };

    }

//  TODO Think about what happends to self rc perfect loops
    std::vector<Sequence> ConstructLoopFromVertex(const KeyWithHash &kh, SequenceBuilder &builder) const {
        DeEdge break_point(kh, origin_.GetUniqueOutgoing(kh));
        Sequence s = ConstructSequenceWithEdge(break_point, builder);
        Kmer kmer = s.start<Kmer>(kmer_size_ + 1) >> 'A';
        for (size_t i = kmer_size_; i < s.size(); i++) {
            kmer = kmer << s[i];
            if (kmer == !kmer)
                return SplitLoop(s, i - kmer_size_);
        }
        return {s};
    }

    void CalculateSequences(kmer_iterator &it,
                            std::vector<Sequence> &sequences) const {
        SequenceBuilder builder;
        std::vector<DeEdge> start_edges;
        start_edges.reserve(8);

        for ( ; it.good(); ++it) {
            KeyWithHash kh = origin_.ConstructKWH(Kmer(kmer_size_, *it));
            AddStartDeEdges(kh, start_edges);

            for (auto edge : start_edges) {
                Sequence s = ConstructSequenceWithEdge(edge, builder);
                if (s < !s)
                    continue;

                sequences.push_back(s);
                TRACE("From " << edge << " calculated sequence\n" << s);
            }
        }
    }

    void CleanCondensed(const Sequence &sequence) {
        Kmer kmer = sequence.start<Kmer>(kmer_size_);
        KeyWithHash kwh = origin_.ConstructKWH(kmer);
        origin_.IsolateVertex(kwh);
        for (size_t pos = kmer_size_; pos < sequence.size(); pos++) {
            kwh = kwh << sequence[pos];
            origin_.IsolateVertex(kwh);
        }
    }

    void CleanCondensed(const std::vector<Sequence> &sequences) {
#       pragma omp parallel for schedule(guided)
        for (size_t i = 0; i < sequences.size(); ++i) {
            CleanCondensed(sequences[i]);
            CleanCondensed(!sequences[i]);
        }
    }

    // This methods collects all loops that were not extracted by finding
    // unbranching paths because there are no junctions on loops.
    const std::vector<Sequence> CollectLoops(unsigned nchunks) {
        INFO("Collecting perfect loops");
        auto its = origin_.kmer_begin(nchunks);
        std::vector<std::vector<KeyWithHash> > starts(its.size());

#       pragma omp parallel for schedule(guided)
        for (size_t i = 0; i < its.size(); ++i) {
            auto &it = its[i];
            for (; it.good(); ++it) {
                KeyWithHash kh = origin_.ConstructKWH(Kmer(kmer_size_, *it));
                if (!IsJunction(kh))
                    starts[i].push_back(kh);
            }
        }

        std::vector<Sequence> result;
        SequenceBuilder builder;
        for (const auto& entry : starts) {
            for (const auto& kwh : entry) {
                if (IsJunction(kwh))
                    continue;

                for (Sequence s : ConstructLoopFromVertex(FindMinimalKMerInLoop(kwh), builder)) {
                    Sequence s_rc = !s;
                    if (s < s_rc)
                        result.push_back(s_rc);
                    else
                        result.push_back(s);

                    CleanCondensed(s);
                    CleanCondensed(s_rc);
                }
            }
        }
        INFO("Collecting perfect loops finished. " << result.size() << " loops collected");
        return result;
    }

public:
    UnbranchingPathExtractor(Index &origin, size_t k)
            : origin_(origin), kmer_size_(k) {}

    //TODO very large vector is returned. But I hate to make all those artificial changes that can fix it.
    const std::vector<Sequence> ExtractUnbranchingPaths(unsigned nchunks) const {
        auto its = origin_.kmer_begin(nchunks);

        INFO("Extracting unbranching paths");
        std::vector<std::vector<Sequence>> sequences(its.size());
#       pragma omp parallel for schedule(guided)
        for (size_t i = 0; i < its.size(); ++i)
            CalculateSequences(its[i], sequences[i]);

        size_t snum = std::accumulate(sequences.begin(), sequences.end(),
                                      0ULL,
                                      [](size_t val, const std::vector<Sequence> &s) {
                                          return val + s.size();
                                      });
        sequences[0].reserve(snum);
        for (size_t i = 1; i < sequences.size(); ++i) {
            sequences[0].insert(sequences[0].end(),
                                std::make_move_iterator(sequences[i].begin()), std::make_move_iterator(sequences[i].end()));
            sequences[i].clear();
            sequences[i].shrink_to_fit();
        }

        INFO("Extracting unbranching paths finished. " << sequences[0].size() << " sequences extracted");
        return sequences[0];
    }

    const std::vector<Sequence> ExtractUnbranchingPathsAndLoops(unsigned nchunks) {
        std::vector<Sequence> result = ExtractUnbranchingPaths(nchunks);
        CleanCondensed(result);
        std::vector<Sequence> loops = CollectLoops(nchunks);
        result.insert(result.end(),
                      std::make_move_iterator(loops.begin()), std::make_move_iterator(loops.end()));
        return result;
    }

private:
    DECL_LOGGER("UnbranchingPathExtractor")
};

template<class Graph>
class FastGraphFromSequencesConstructor {
private:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef RtSeq Kmer;
    typedef kmers::DeBruijnExtensionIndex<> Index;
    size_t kmer_size_;
    Index &origin_;

    class LinkRecord {
    private:
        uint64_t hash_and_mask_;
        EdgeId edge_;

        size_t BitBool(bool flag) const {
            if (flag)
                return 1;
            return 0;
        }

    public:
        uint64_t GetHash() const { return hash_and_mask_ >> 2; }
        bool IsRC() const { return hash_and_mask_ & 2; }
        bool IsStart() const { return hash_and_mask_ & 1; }
        EdgeId GetEdge() const { return edge_; }
        bool IsInvalid() const { return hash_and_mask_ + 1 == 0 && edge_ == EdgeId(); }

        LinkRecord(uint64_t hash, EdgeId edge, bool is_start, bool is_rc)
                : hash_and_mask_((hash << 2) | (BitBool(is_rc) << 1)| BitBool(is_start)), edge_(edge) { }

        LinkRecord()
                : hash_and_mask_(-1ul) {}

        uint64_t EdgeAndMask() const {
            return (edge_.int_id() << 2) | (BitBool(IsRC()) << 1) | BitBool(IsStart());
        }

        static bool CompareByVertexKMerEdgeIdAndMask(const LinkRecord &r1, const LinkRecord& r2) {
            return std::make_tuple(r1.GetHash(), r1.EdgeAndMask()) < std::make_tuple(r2.GetHash(), r2.EdgeAndMask());
        }
    };

    LinkRecord StartLink(const EdgeId &edge, const Sequence &sequence) const {
        Kmer kmer(kmer_size_, sequence);
        Kmer kmer_rc = !kmer;
        if (kmer < kmer_rc)
            return LinkRecord(origin_.ConstructKWH(kmer).idx(), edge, true, false);
        else
            return LinkRecord(origin_.ConstructKWH(kmer_rc).idx(), edge, true, true);
    }

    LinkRecord EndLink(const EdgeId &edge, const Sequence &sequence) const {
        Kmer kmer(kmer_size_, sequence, sequence.size() - kmer_size_);
        Kmer kmer_rc = !kmer;
        if (kmer < kmer_rc)
            return LinkRecord(origin_.ConstructKWH(kmer).idx(), edge, false, false);
        else
            return LinkRecord(origin_.ConstructKWH(kmer_rc).idx(), edge, false, true);
    }

    void CollectLinkRecords(typename Graph::HelperT &helper, const Graph &graph,
                            std::vector<LinkRecord> &records, const std::vector<Sequence> &sequences) const {
        size_t size = sequences.size();
        uint64_t min_id = graph.min_id();
        records.resize(size * 2, LinkRecord(0, EdgeId(), false, false));
#       pragma omp parallel for schedule(guided)
        for (size_t i = 0; i < size; ++i) {
            size_t j = i << 1;
            EdgeId edge = helper.AddEdge(DeBruijnEdgeData(sequences[i]), min_id + j);
            records[j] = StartLink(edge, sequences[i]);
            if (graph.conjugate(edge) != edge)
                records[j + 1] = EndLink(edge, sequences[i]);
            else
                records[j + 1] = LinkRecord();
        }
    }

    void LinkEdge(typename Graph::HelperT &helper, const Graph &graph, const VertexId v,
                  const EdgeId edge, const bool is_start, const bool is_rc) const {
        VertexId v1 = v;
        if (is_rc)
            v1 = graph.conjugate(v);

        if (is_start)
            helper.LinkOutgoingEdge(v1, edge);
        else
            helper.LinkIncomingEdge(v1, edge);
    }

public:
    FastGraphFromSequencesConstructor(size_t k, Index &origin)
            : kmer_size_(k), origin_(origin) {}

    void ConstructGraph(Graph &graph, const std::vector<Sequence> &sequences) const {
        typename Graph::HelperT helper = graph.GetConstructionHelper();

        std::vector<LinkRecord> records;
        size_t seq_size = sequences.size();
        INFO("Total " << 2*seq_size << " edges to create");
        graph.ereserve(2*seq_size + seq_size / 100);

        INFO("Collecting link records")
        CollectLinkRecords(helper, graph, records, sequences);
        INFO("Ordering link records")
        // We sort by Vertex and then by EdgeID and RC/Start mask in order to combine together records accociated with the same vertex with a special order in each group
        parallel::sort(records.begin(), records.end(), LinkRecord::CompareByVertexKMerEdgeIdAndMask);
        INFO("Sorting done");

        // Now we extract starting positions of each vertex group
        std::vector<size_t> unique_record_indices;
        for (size_t i = 0; i < records.size(); i++) {
            if (i == 0 || records[i].GetHash() != records[i - 1].GetHash()) {
                if (!records[i].IsInvalid()) {
                    unique_record_indices.push_back(i);
                }
            }
        }
        // Now we sort vertices by their lowest edge and mask (they are unique since each edge has only one start and one stop).
        // It is a deterministic order while ordering by vertex kmer perfect hash is not (hashes are dependent on nthreads/nnodes)
        INFO("Sorting LinkRecords...");
        parallel::sort(unique_record_indices.begin(), unique_record_indices.end(),
                       [&records](size_t i, size_t j) { return records[i].EdgeAndMask() < records[j].EdgeAndMask(); });
        INFO("LinkRecords sorted");
        size_t size = unique_record_indices.size();
        INFO("Total " << size << " vertices to create");
        graph.vreserve(2 * size + size / 100);

        INFO("Connecting the graph");
        uint64_t min_id = graph.min_id();
#       pragma omp parallel for schedule(guided)
        for (size_t vertex_num = 0; vertex_num < size; ++vertex_num) {
            size_t i = unique_record_indices[vertex_num];

            VertexId v = helper.CreateVertex(DeBruijnVertexData(graph.k()), min_id + (vertex_num << 1));
            for (size_t j = i; j < records.size() && records[j].GetHash() == records[i].GetHash(); j++) {
                LinkEdge(helper, graph, v, records[j].GetEdge(), records[j].IsStart(), records[j].IsRC());
            }
        }
    }
};

/*
 * Constructs DeBruijnGraph from DeBruijnExtensionIndex using "new DeBruijnGraphExtentionConstructor(DeBruijn).ConstructGraph(DeBruijnGraph, Index)"
 */
template<class Graph>
class DeBruijnGraphExtentionConstructor {
private:
    typedef typename Graph::EdgeId EdgeId;
    typedef kmers::DeBruijnExtensionIndex<> DeBruijn;
    typedef typename Graph::VertexId VertexId;
    typedef RtSeq Kmer;

    Graph &graph_;
    DeBruijn &origin_;
    size_t kmer_size_;

public:
    DeBruijnGraphExtentionConstructor(Graph& graph, DeBruijn &origin) :
            graph_(graph), origin_(origin), kmer_size_(graph.k()) {
    }

    void ConstructGraph(bool keep_perfect_loops) {
        std::vector<Sequence> edge_sequences;
        unsigned nchunks = 16 * omp_get_max_threads();
        if (keep_perfect_loops)
            edge_sequences = UnbranchingPathExtractor(origin_, kmer_size_).ExtractUnbranchingPathsAndLoops(nchunks);
        else
            edge_sequences = UnbranchingPathExtractor(origin_, kmer_size_).ExtractUnbranchingPaths(nchunks);
        INFO("Sorting edges...");
        parallel::sort(edge_sequences.begin(), edge_sequences.end(), Sequence::RawCompare);
        INFO("Edges sorted");
        FastGraphFromSequencesConstructor<Graph>(kmer_size_, origin_).ConstructGraph(graph_, edge_sequences);
    }

private:
    DECL_LOGGER("DeBruijnGraphConstructor")
};

}
