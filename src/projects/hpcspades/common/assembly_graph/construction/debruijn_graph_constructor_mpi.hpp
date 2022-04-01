#pragma once
//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "projects/hpcspades/common/pipeline/partask_mpi.hpp"
#include "io/binary/graph.hpp"
#include "common/assembly_graph/construction/debruijn_graph_constructor.hpp"

namespace debruijn_graph {
template<class Graph>
class DeBruijnGraphExtentionConstructorTask {
 private:
    typedef typename Graph::EdgeId EdgeId;
    typedef kmers::DeBruijnExtensionIndex<> Index;
    typedef typename Graph::VertexId VertexId;
    typedef RtSeq Kmer;

    bool collect_loops_;

 public:
    DeBruijnGraphExtentionConstructorTask(std::istream &is) {
        io::binary::BinRead(is, collect_loops_);
    }

    DeBruijnGraphExtentionConstructorTask(bool collect_loops) : collect_loops_{collect_loops} {}

    std::ostream &serialize(std::ostream &os) const {
        io::binary::BinWrite(os, collect_loops_);
        return os;
    }

    template<typename... Args>
    auto make_splitter(size_t size, Args &&...) {
        return partask::make_seq_plus_n_generator(size);
    }

    void process(std::istream &is, std::ostream &os, Graph &g, Index &index) {
        size_t n = 0;
        std::vector<size_t> chunks = partask::get_seq_plus_n(is, n);
        if (!chunks.size()) {
            INFO("Empty job, skipping");
        }

        auto iters = index.kmer_begin(n);

        std::vector<typename Index::kmer_iterator> local_iters;
        for (size_t i : chunks) {
            if (i < iters.size()) {
                local_iters.push_back(std::move(iters[i]));
            }
        }

        UnbranchingPathExtractor extractor(index, g.k());
        auto seqs = extractor.ExtractUnbranchingPaths(local_iters);
        io::binary::BinWrite(os, partask::fast_local_transfer(seqs));
    }

    void merge(const std::vector<std::istream *> &piss, Graph &g, Index &index) {
        std::vector<Sequence> seqs;
        for (size_t i = 0; i < piss.size(); ++i) {
            auto &is = *piss[i];
            if (is.peek() != EOF) {
                std::vector<Sequence> local_seqs;
                io::binary::BinRead(is, partask::fast_local_transfer(local_seqs));
                seqs.insert(seqs.end(),
                            std::make_move_iterator(local_seqs.begin()), std::make_move_iterator(local_seqs.end()));
            }
        }

        {
            TIME_TRACE_SCOPE("RemoveSequences");
            index.RemoveSequences(seqs);
        }

        if (collect_loops_) {
            UnbranchingPathExtractor extractor(index, g.k());
            std::vector<Sequence> loops = extractor.CollectLoops(omp_get_max_threads());
            seqs.insert(seqs.end(),
                        std::make_move_iterator(loops.begin()), std::make_move_iterator(loops.end()));
        }

        INFO("Sorting edges...");
        {
            TIME_TRACE_SCOPE("Sorting edges");
            parallel::sort(seqs.begin(), seqs.end(), Sequence::RawCompare);
        }
        INFO("Sorting edges finished");

        FastGraphFromSequencesConstructor<Graph>(g.k(), index).ConstructGraph(g, seqs);
    }

    void sync(Graph &g, Index &) {
        auto serialize = [](std::ostream &os, const Graph &g) {
            io::binary::GraphIO<Graph>().BinWrite(os, g);
        };
        auto deserialize = [](std::istream &is, Graph &g) {
            io::binary::GraphIO<Graph>().BinRead(is, g);
        };
        partask::broadcast(g, serialize, deserialize);
    }
};
} //namespace debruijn_graph
