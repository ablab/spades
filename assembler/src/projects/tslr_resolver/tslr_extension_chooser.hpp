#pragma once

#include <algorithms/path_extend/extension_chooser.hpp>
#include <algorithms/path_extend/path_filter.hpp>
#include <algorithms/path_extend/path_extender.hpp>
#include <algorithms/path_extend/pe_resolver.hpp>
#include <algorithms/path_extend/path_extend_launch.hpp>
#include <assembly_graph/graph_support/scaff_supplementary.hpp>
#include <cmath>
#include "barcode_mapper.hpp"

using namespace path_extend;

namespace tslr_resolver {

    template <class Graph> //TODO remove from here
        bool IsEdgeUnique(const Graph& g, const EdgeId& edge) {
            return g.coverage(edge) < 1.1 * cfg::get().ts_res.reference_cov;
        }

        class TrivialTSLRExtensionChooser : public ExtensionChooser {

            shared_ptr<BarcodeMapper> bmapper_;
            const EdgesPositionHandler<Graph>& edge_pos_;
            int reference_cov_;
            size_t len_threshold_;
            double relative_diff_threshold_;

        public:
            TrivialTSLRExtensionChooser(const conj_graph_pack& gp, const int& reference_cov, const size_t& len_threshold,
                const double& relative_diff_threshold) :
                    ExtensionChooser(gp.g), bmapper_(gp.barcode_mapper), 
                    edge_pos_(gp.edge_pos), reference_cov_(reference_cov), 
                    len_threshold_(len_threshold), relative_diff_threshold_(relative_diff_threshold) {
            }

            EdgeContainer Filter(const BidirectionalPath &path, const EdgeContainer &edges) const override {
                auto result = EdgeContainer();
                if (edges.size() == 0) {
                    return result;
                }
                //Find long unique edge earlier in path
                //FIXME code duplication
                bool long_single_edge_exists = false;
                EdgeId decisive_edge;
                for (int i = static_cast<int> (path.Size()) - 1; !long_single_edge_exists && i >= 0; --i) {
                    EdgeId current_edge = path[i];
                    if (IsEdgeUnique<Graph>(g_, current_edge) &&
                            g_.length(current_edge) > len_threshold_) {
                        long_single_edge_exists = true;
                        decisive_edge = current_edge;
                    }
                }

                //Exclude decisive edge from candidates
                auto edges_copy = edges;
                EraseEdge(edges_copy, decisive_edge);
                if (edges_copy.size() == 1) {
                    return edges_copy;
                }
                if (!long_single_edge_exists || edges_copy.size() == 0) {
                    if (edges_copy.size() == 0) {
                        DEBUG("Only decisive edge found");
                    }
                    return result;
                }
                
                //Find edges with best barcode score
                auto fittest_edge = *(std::max_element(edges_copy.begin(), edges_copy.end(),
                                                     [this, & decisive_edge](const EdgeWithDistance& edge1, const EdgeWithDistance& edge2) {
                                                         return this->bmapper_->IntersectionSizeNormalizedBySecond(decisive_edge, edge1.e_) <
                                                                this->bmapper_->IntersectionSizeNormalizedBySecond(decisive_edge, edge2.e_);
                                                     }));
                double best_score = bmapper_->IntersectionSizeNormalizedBySecond(decisive_edge, fittest_edge.e_);

                DEBUG("fittest edge " << fittest_edge.e_.int_id());
                DEBUG("score " << best_score);

                std::vector <EdgeWithDistance> best_candidates;
                std::copy_if(edges_copy.begin(), edges_copy.end(), std::back_inserter(best_candidates), 
                                    [this, &decisive_edge, best_score](const EdgeWithDistance& edge) {
                                        return this->bmapper_->IntersectionSizeNormalizedBySecond(decisive_edge, edge.e_) +
                                        relative_diff_threshold_ > best_score &&
                                        IsEdgeUnique<Graph>(g_, edge.e_);
                                    });

                if (best_candidates.size() == 1) {
                    result.push_back(fittest_edge);
                    return result;
                }
                if (best_candidates.size() == 0) {
                    return result;
                }
                //Try to find topologically closest edge
                DEBUG("Several candidates found. Further filtering.");
                std::nth_element(edges_copy.begin(), edges_copy.begin() + 1, edges_copy.end(), 
                                                 [this, & decisive_edge](const EdgeWithDistance& edge1, const EdgeWithDistance& edge2) {
                                                     return this->bmapper_->IntersectionSizeNormalizedBySecond(decisive_edge, edge1.e_) >
                                                            this->bmapper_->IntersectionSizeNormalizedBySecond(decisive_edge, edge2.e_);
                                                 });
                auto second_best_edge = edges_copy.at(1);
                double second_best_score = bmapper_->IntersectionSizeNormalizedBySecond(decisive_edge, second_best_edge.e_);
                DEBUG("Second best edge " << second_best_edge.e_.int_id());
                DEBUG("second best score " << second_best_score);
                DEBUG(best_candidates.size() << " best candidates");

                auto it = FindClosestEdge(best_candidates);
                if (it == best_candidates.end()) {
                    DEBUG("Closest edge wasn't found");
                    for (auto edge : best_candidates) {
                        result.push_back(edge);
                    }
                }
                else {
                    DEBUG("Found topologically closest edge");
                    result.push_back(*it);
                }
                return result;
            }

        private:

            void EraseEdge (EdgeContainer& edges, const EdgeId& edge_to_be_erased) const {
                size_t ind = edges.size() + 1;
                for (size_t i = 0; i < edges.size(); ++i) {
                    if (edges[i].e_ == edge_to_be_erased) {
                        ind = i;
                        break;
                    }
                }
                if (ind != edges.size() + 1)
                    edges.erase(edges.begin() + ind);
            }

            vector<EdgeWithDistance>::const_iterator FindClosestEdge(const vector<EdgeWithDistance>& edges) const {
                //Make it more effective if needed
                auto it = edges.begin();
                do {
                    auto edge = *it;
                    size_t path_len_bound = cfg::get().ts_res.topsort_bound;
                    VertexId start_vertex = g_.EdgeEnd(edge.e_);
                    auto dijkstra = DijkstraHelper<Graph>::CreateBoundedDijkstra(g_, path_len_bound);
                    dijkstra.Run(start_vertex);
                    bool can_reach_everyone = true;
                    for (auto other_edge: edges) {
                        if (other_edge.e_ == edge.e_) 
                            continue;
                        auto other_start = g_.EdgeStart(other_edge.e_);
                        if (!dijkstra.DistanceCounted(other_start)) {
                            can_reach_everyone = false;
                            break;
                        }
                    }
                    if (can_reach_everyone) {
                        return it;
                    }
                    ++it;
                } while (it != edges.end());
                return it;
            }
            DECL_LOGGER("TslrExtensionChooser")
        };
} //tslr_resolver