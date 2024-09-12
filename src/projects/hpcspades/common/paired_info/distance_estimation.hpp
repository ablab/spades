//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef MPI_DISTANCE_ESTIMATION_HPP_
#define MPI_DISTANCE_ESTIMATION_HPP_

#include "common/paired_info/distance_estimation.hpp"
#include "projects/hpcspades/common/pipeline/partask_mpi.hpp"

namespace omnigraph {
    namespace de {
        class DistanceEstimatorMPI : public DistanceEstimator {
            typedef DistanceEstimator base;
            typedef std::vector <size_t> GraphLengths;
            typedef std::vector <std::pair<int, double>> EstimHist;
            typedef std::pair <debruijn_graph::EdgeId, debruijn_graph::EdgeId> EdgePair;

        protected:
            typedef typename base::InPairedIndex InPairedIndex;
            typedef typename base::OutPairedIndex OutPairedIndex;
            typedef typename base::InHistogram InHistogram;
            typedef typename base::OutHistogram OutHistogram;

        public:
            DistanceEstimatorMPI(const debruijn_graph::Graph &graph,
                                 const InPairedIndex &index,
                                 const GraphDistanceFinder &distance_finder,
                                 const PairInfoChecker &checker,
                                 size_t linkage_distance, size_t max_distance,
                                 std::unique_ptr<DistanceEstimator> base_dist_estimator)
                    : base(graph, index, distance_finder, checker, linkage_distance, max_distance),
                      dist_estimator_(std::move(base_dist_estimator)) {}

            class DistanceEstimatorTask {
                DistanceEstimatorTask() = default;

            public:
                DistanceEstimatorTask(std::vector <debruijn_graph::EdgeId> &edges,
                                      unsigned int nthreads) : edges_(edges), nthreads_(nthreads) {};

                DistanceEstimatorTask(std::istream &is) {
                    io::binary::BinRead(is, edges_, nthreads_);

                }

                std::ostream &serialize(std::ostream &os) const {
                    io::binary::BinWrite(os, edges_, nthreads_);
                    return os;
                }

                auto make_splitter(size_t, const InPairedIndex &, const DistanceEstimator &,
                                   PairedInfoIndexT <debruijn_graph::Graph> & /*result*/) {
                    return partask::make_seq_along_generator(edges_);
                }

                void process(std::istream &is, std::ostream &os, const InPairedIndex &index,
                             const DistanceEstimator &self, PairedInfoIndexT <debruijn_graph::Graph> & /*result*/) {
                    DEBUG("Processing");
                    auto edges_id = partask::get_seq(is);
                    Buffer buffer(self.graph());
#                   pragma omp parallel for num_threads(nthreads_) schedule(guided, 10)
                    for (size_t i = 0; i < edges_id.size(); ++i) {
                        debruijn_graph::EdgeId edge = edges_[edges_id[i]];
                        self.ProcessEdge(edge, index, buffer);
                    }

                    buffer.BinWrite(os);
                }

                auto merge(const std::vector<std::istream *> &piss,
                           const InPairedIndex &,
                           const DistanceEstimator &self,
                           PairedInfoIndexT <debruijn_graph::Graph> &result) {
                    for (auto pis: piss) {
                        Buffer buffer(self.graph());
                        buffer.BinRead(*pis);
                        result.MergeAssign(buffer);
                    }
                }

            private:
                std::vector <debruijn_graph::EdgeId> edges_;
                unsigned nthreads_;
            };

            void Estimate(OutPairedIndex &result, size_t nthreads) const override;

            friend DistanceEstimatorTask;
        private:
            std::unique_ptr<DistanceEstimator> dist_estimator_;

            const std::string Name() const override {
                const std::string my_name = dist_estimator_->Name() + "_MPI";
                return my_name;
            }

            DECL_LOGGER("DistanceEstimatorMPI");
        };
    }
}

#endif /* MPI_DISTANCE_ESTIMATION_HPP_ */
