#pragma once
//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/core/construction_helper.hpp"
#include "utils/standard_base.hpp"
#include "utils/indices/kmer_extension_index.hpp"
#include "utils/openmp_wrapper.h"
#include "utils/parallel_wrapper.hpp"

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

#   pragma omp parallel for schedule(guided)
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
            queueSize = min(size_t(double(queueSize) * queueGrowthRate), queueMaxSize);
        }
    }

private:
    DECL_LOGGER("DeBruijnGraphConstructor")
};

class UnbranchingPathFinder {
private:
    typedef DeBruijnExtensionIndex<> Index;
    typedef RtSeq Kmer;
    typedef Index::kmer_iterator kmer_iterator;
    typedef Index::KeyWithHash KeyWithHash;
    typedef Index::DeEdge DeEdge;

    Index &origin_;
    size_t kmer_size_;

public:
    UnbranchingPathFinder(Index &origin, size_t kmer_size) : origin_(origin), kmer_size_(kmer_size) {
    }

    bool StepRightIfPossible(DeEdge &edge) {
        if (origin_.CheckUniqueOutgoing(edge.end) && origin_.CheckUniqueIncoming(edge.end)) {
            edge = DeEdge(edge.end, origin_.GetUniqueOutgoing(edge.end));
            return true;
        }
        return false;
    }

    Sequence ConstructSeqGoingRight(DeEdge edge) {
        SequenceBuilder s;
        s.append(edge.start.key());
        s.append(edge.end[kmer_size_ - 1]);
        DeEdge initial = edge;
        while (StepRightIfPossible(edge) && edge != initial) {
            s.append(edge.end[kmer_size_ - 1]);
        }
        return s.BuildSequence();
    }

    Sequence ConstructSequenceWithEdge(DeEdge edge) {
        return ConstructSeqGoingRight(edge);
    }

    //Loop consists of 4 parts: 2 selfRC k+1-mers and two sequences of arbitrary length RC to each other; pos is a position of one of selfRC edges
    vector<Sequence> SplitLoop(Sequence s, size_t pos) {
        return {s.Subseq(pos, pos + kmer_size_ + 1), s.Subseq(pos + 1, s.size() - kmer_size_) + s.Subseq(0, pos + kmer_size_)};

    }

//TODO Think about what happends to self rc perfect loops
    vector<Sequence> ConstructLoopFromVertex(const KeyWithHash &kh) {
        DeEdge break_point(kh, origin_.GetUniqueOutgoing(kh));
        Sequence s = ConstructSequenceWithEdge(break_point);
        Kmer kmer = s.start<Kmer>(kmer_size_ + 1) >> 'A';
        for(size_t i = kmer_size_; i < s.size(); i++) {
            kmer = kmer << s[i];
            if (kmer == !kmer) {
                return SplitLoop(s, i - kmer_size_);
            }
        }
        return {s};
    }
};

class UnbranchingPathExtractor {
private:
    typedef DeBruijnExtensionIndex<> Index;
    typedef RtSeq Kmer;
    typedef Index::kmer_iterator kmer_iterator;
    typedef Index::DeEdge DeEdge;
    typedef Index::KeyWithHash KeyWithHash;

    Index &origin_;
    size_t kmer_size_;

    bool IsJunction(KeyWithHash kh) const {
        return !(origin_.CheckUniqueOutgoing(kh) && origin_.CheckUniqueIncoming(kh));
    }

    void AddStartDeEdgesForVertex(KeyWithHash kh, std::vector<DeEdge>& start_edges) const {
        for (char next = 0; next < 4; next++) {
            if (origin_.CheckOutgoing(kh, next)) {
                TRACE("Added to queue " << DeEdge(kh, origin_.GetOutgoing(kh, next)));
                start_edges.push_back(DeEdge(kh, origin_.GetOutgoing(kh, next)));
            }
        }
    }

    void AddStartDeEdges(kmer_iterator &it, size_t queueSize,
                  std::vector<DeEdge>& start_edges) const {
        for (; start_edges.size() < queueSize && it.good(); ++it) {
            KeyWithHash kh = origin_.ConstructKWH(Kmer(kmer_size_, *it));
            if (IsJunction(kh)) {
                AddStartDeEdgesForVertex(kh, start_edges);
                KeyWithHash kh_inv = !kh;
                if(!(kh_inv.is_minimal())) {
                    AddStartDeEdgesForVertex(kh_inv, start_edges);
                }
            }
        }
    }

    void CalculateSequences(std::vector<DeEdge> &edges,
                            std::vector<Sequence> &sequences, UnbranchingPathFinder &finder) const {
        size_t size = edges.size();
        size_t start = sequences.size();
        sequences.resize(start + size);

#   pragma omp parallel for schedule(guided)
        for (size_t i = 0; i < size; ++i) {
            sequences[start + i] = finder.ConstructSequenceWithEdge(edges[i]);
            TRACE("From " << edges[i] << " calculated sequence");
            TRACE(sequences[start + i]);
        }
    }

    void CleanCondensed(const Sequence &sequence) {
        Kmer kmer = sequence.start<Kmer>(kmer_size_);
        KeyWithHash kwh = origin_.ConstructKWH(kmer);
        origin_.IsolateVertex(kwh);
        for(size_t pos = kmer_size_; pos < sequence.size(); pos++) {
            kwh = kwh << sequence[pos];
            origin_.IsolateVertex(kwh);
        }
    }

    void CleanCondensed(const std::vector<Sequence> &sequences) {
#   pragma omp parallel for schedule(guided)
        for (size_t i = 0; i < sequences.size(); ++i) {
            CleanCondensed(sequences[i]);
        }
    }

    //This methods collects all loops that were not extracted by finding unbranching paths because there are no junctions on loops.
    //TODO make parallel
    const std::vector<Sequence> CollectLoops() {
        INFO("Collecting perfect loops");
        UnbranchingPathFinder finder(origin_, kmer_size_);
        std::vector<Sequence> result;
        for (kmer_iterator it = origin_.kmer_begin(); it.good(); ++it) {
            KeyWithHash kh = origin_.ConstructKWH(Kmer(kmer_size_, *it));
            if (!IsJunction(kh)) {
                vector<Sequence> loop = finder.ConstructLoopFromVertex(kh);
                for(Sequence s: loop) {
                    result.push_back(s);
                    CleanCondensed(s);
                    if(s != (!s)) {
                        result.push_back(!s);
                    }
                }
            }
        }
        INFO("Collecting perfect loops finished. " << result.size() << " loops collected");
        return result;
    }

public:
    UnbranchingPathExtractor(Index &origin, size_t k) : origin_(origin), kmer_size_(k) {
    }

    //TODO very large vector is returned. But I hate to make all those artificial changes that can fix it.
    const std::vector<Sequence> ExtractUnbranchingPaths(size_t queueMinSize, size_t queueMaxSize,
                                                        double queueGrowthRate) {
        INFO("Extracting unbranching paths");
        UnbranchingPathFinder finder(origin_, kmer_size_);
        std::vector<Sequence> result;
        size_t queueSize = queueMinSize;
        std::vector<DeEdge> start_edges;
        std::vector<Sequence> sequences;
        start_edges.reserve(queueSize);
        auto it = origin_.kmer_begin();
        while (it.good()) {
            AddStartDeEdges(it, queueSize, start_edges); // format a queue of junction kmers
            CalculateSequences(start_edges, sequences, finder); // in parallel
            start_edges.clear();
            queueSize = min((size_t) ((double) queueSize * queueGrowthRate), queueMaxSize);
        }
        INFO("Extracting unbranching paths finished. " << sequences.size() << " sequences extracted");
        return sequences;
    }

    const std::vector<Sequence> ExtractUnbranchingPathsAndLoops(size_t queueMinSize, size_t queueMaxSize,
                                                                double queueGrowthRate) {
        std::vector<Sequence> result = ExtractUnbranchingPaths(queueMinSize, queueMaxSize, queueGrowthRate);
        CleanCondensed(result);
        std::vector<Sequence> loops = CollectLoops();
        for(auto it = loops.begin(); it != loops.end(); ++it) {
            result.push_back(*it);
        }
        return result;
    }

private:
    DECL_LOGGER("UnbranchingPathExtractor")
};

/*
 * Only works for Conjugate dbg
 */
template<class Graph>
class FastGraphFromSequencesConstructor {
private:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef RtSeq Kmer;
    typedef DeBruijnExtensionIndex<> Index;
    size_t kmer_size_;
    Index &origin_;

    class LinkRecord {
    private:
        size_t hash_and_mask_;
        EdgeId edge_;

        size_t BitBool(bool flag) const {
            if(flag)
                return 1;
            return 0;
        }

    public:
        size_t GetHash() const {
            return hash_and_mask_ >> 2;
        }

        bool IsRC() const {
            return hash_and_mask_ & 2;
        }

        bool IsStart() const {
            return hash_and_mask_ & 1;
        }


        EdgeId GetEdge() const {
            return edge_;
        }

        LinkRecord(size_t hash, EdgeId edge, bool is_start, bool is_rc) :
                hash_and_mask_((hash << 2) | (BitBool(is_rc) << 1)| BitBool(is_start)), edge_(edge) {
        }

        LinkRecord() :
                hash_and_mask_(-1ul), edge_(0) {
        }

        bool IsInvalid() {
            return hash_and_mask_ + 1 == 0 && edge_ == EdgeId(0);
        }

        bool operator<(const LinkRecord &other) const {
            if(this->hash_and_mask_ == other.hash_and_mask_)
                return this->edge_ < other.edge_;
            return this->hash_and_mask_ < other.hash_and_mask_;
        }
    };

    LinkRecord StartLink(const EdgeId &edge, const Sequence &sequence) const {
        Kmer kmer(kmer_size_, sequence);
        Kmer kmer_rc = !kmer;
        if(kmer < kmer_rc)
            return LinkRecord(origin_.ConstructKWH(kmer).idx(), edge, true, false);
        else
            return LinkRecord(origin_.ConstructKWH(kmer_rc).idx(), edge, true, true);
    }

    LinkRecord EndLink(const EdgeId &edge, const Sequence &sequence) const {
        Kmer kmer(kmer_size_, sequence, sequence.size() - kmer_size_);
        Kmer kmer_rc = !kmer;
        if(kmer < kmer_rc)
            return LinkRecord(origin_.ConstructKWH(kmer).idx(), edge, false, false);
        else
            return LinkRecord(origin_.ConstructKWH(kmer_rc).idx(), edge, false, true);
    }

    void CollectLinkRecords(typename Graph::HelperT &helper, const Graph &graph, vector<LinkRecord> &records, const vector<Sequence> &sequences) const {
        size_t size = sequences.size();
        records.resize(size * 2, LinkRecord(0, EdgeId(0), false, false));
        restricted::IdSegmentStorage id_storage = helper.graph().GetGraphIdDistributor().Reserve(size * 2);
#   pragma omp parallel for schedule(guided)
        for (size_t i = 0; i < size; ++i) {
            size_t j = i << 1;
            auto id_distributor = id_storage.GetSegmentIdDistributor(j, j + 2);//indices for two edges are required
            EdgeId edge = helper.AddEdge(DeBruijnEdgeData(sequences[i]), id_distributor);
            records[j] = StartLink(edge, sequences[i]);
            if(graph.conjugate(edge) != edge)
                records[j + 1] = EndLink(edge, sequences[i]);
            else
                records[j + 1] = LinkRecord();
        }
    }

    void LinkEdge(typename Graph::HelperT &helper, const Graph &graph, const VertexId v, const EdgeId edge, const bool is_start, const bool is_rc) const {
        VertexId v1 = v;
        if(is_rc) {
            v1 = graph.conjugate(v);
        }
        if(is_start) {
            helper.LinkOutgoingEdge(v1, edge);
        } else {
            helper.LinkIncomingEdge(v1, edge);
        }
    }

public:
    FastGraphFromSequencesConstructor(size_t k, Index &origin) : kmer_size_(k), origin_(origin) {
    }

    void ConstructGraph(Graph &graph, const vector<Sequence> &sequences) const {
        typename Graph::HelperT helper = graph.GetConstructionHelper();
        vector<LinkRecord> records;
        CollectLinkRecords(helper, graph, records, sequences);//TODO make parallel
        parallel::sort(records.begin(), records.end());
        size_t size = records.size();
        vector<vector<VertexId>> vertices_list(omp_get_max_threads());
        restricted::IdSegmentStorage id_storage = helper.graph().GetGraphIdDistributor().Reserve(size * 2);
#   pragma omp parallel for schedule(guided)
        for(size_t i = 0; i < size; i++) {
            if(i != 0 && records[i].GetHash() == records[i - 1].GetHash()) {
                continue;
            }
            if(records[i].IsInvalid())
                continue;
            auto id_distributor = id_storage.GetSegmentIdDistributor(i << 1, (i << 1) + 2);
            VertexId v = helper.CreateVertex(DeBruijnVertexData(), id_distributor);
            vertices_list[omp_get_thread_num()].push_back(v);
            for(size_t j = i; j < size && records[j].GetHash() == records[i].GetHash(); j++) {
                LinkEdge(helper, graph, v, records[j].GetEdge(), records[j].IsStart(), records[j].IsRC());
            }
        }
        for(size_t i = 0; i < vertices_list.size(); i++)
            helper.AddVerticesToGraph(vertices_list[i].begin(), vertices_list[i].end());
    }
};

/*
 * Constructs DeBruijnGraph from DeBruijnExtensionIndex using "new DeBruijnGraphExtentionConstructor(DeBruijn).ConstructGraph(DeBruijnGraph, Index)"
 */
template<class Graph>
class DeBruijnGraphExtentionConstructor {
private:
    typedef typename Graph::EdgeId EdgeId;
    typedef DeBruijnExtensionIndex<> DeBruijn;
    typedef typename Graph::VertexId VertexId;
    typedef RtSeq Kmer;

    Graph &graph_;
    DeBruijn &origin_;
    size_t kmer_size_;

    void FilterRC(std::vector<Sequence> &edge_sequences) {
        size_t size = 0;
        for(size_t i = 0; i < edge_sequences.size(); i++) {
            if(!(edge_sequences[i] < !edge_sequences[i])) {
                edge_sequences[size] = edge_sequences[i];
                size++;
            }
        }
        edge_sequences.resize(size);
    }

public:
    DeBruijnGraphExtentionConstructor(Graph& graph, DeBruijn &origin) :
            graph_(graph), origin_(origin), kmer_size_(graph.k()) {
    }

    void ConstructGraph(size_t queueMinSize, size_t queueMaxSize,
            double queueGrowthRate, bool keep_perfect_loops) {
        std::vector<Sequence> edge_sequences;
        if(keep_perfect_loops)
            edge_sequences = UnbranchingPathExtractor(origin_, kmer_size_).ExtractUnbranchingPathsAndLoops(queueMinSize, queueMaxSize, queueGrowthRate);
        else
            edge_sequences = UnbranchingPathExtractor(origin_, kmer_size_).ExtractUnbranchingPaths(queueMinSize, queueMaxSize, queueGrowthRate);
        FilterRC(edge_sequences);
        FastGraphFromSequencesConstructor<Graph>(kmer_size_, origin_).ConstructGraph(graph_, edge_sequences);
    }

private:
    DECL_LOGGER("DeBruijnGraphConstructor")
};

}
