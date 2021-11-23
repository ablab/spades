//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/stage.hpp"
#include "pipeline/mpi_stage.hpp"
#include "assembly_graph/core/graph.hpp"
#include "alignment/sequence_mapper_notifier.hpp"

namespace debruijn_graph {
    namespace mismatches {
        struct MismatchEdgeInfo;
        class MismatchStatistics;

        class MismatchShallNotPass {
        private:
            typedef typename Graph::EdgeId EdgeId;
            typedef typename Graph::VertexId VertexId;
            typedef std::function<void(SequenceMapperListener *, const SequenceMapper<Graph> &,
                                       io::ReadStreamList<io::SingleReadSeq> &streams)> ProccessLibFuncT;

            graph_pack::GraphPack &gp_;
            Graph &graph_;
            const size_t k_;
            const double relative_threshold_;
            const ProccessLibFuncT &proccess_lib_func_;
            const size_t num_readers_;

            EdgeId CorrectNucl(EdgeId edge, size_t position, char nucl);

            EdgeId CorrectNucls(EdgeId edge, const std::vector<std::pair<size_t, char>> &mismatches);

            std::vector<std::pair<size_t, char>> FindMismatches(EdgeId edge, const MismatchEdgeInfo &statistics) const;

            size_t CorrectEdge(EdgeId edge, const MismatchEdgeInfo &statistics);

            size_t CorrectAllEdges(const MismatchStatistics &statistics);

            size_t ParallelStopMismatchIteration();

        public:
            MismatchShallNotPass(const ProccessLibFuncT &processLib, graph_pack::GraphPack &gp,
                                 double relative_threshold = 1.5,
                                 size_t num_readers = 0);


            size_t ParallelStopAllMismatches(size_t max_iterations);
        };
    }


    class MismatchCorrection : public spades::MPIAssemblyStage {
    public:
        MismatchCorrection()
                : MPIAssemblyStage("Mismatch Correction", "mismatch_correction") {}

        void run(graph_pack::GraphPack &gp, const char *) override;
    };
}

