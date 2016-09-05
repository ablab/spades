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
        class TrivialTSLRExtensionChooser : public ExtensionChooser {

            shared_ptr<BarcodeMapper> bmapper_;
            size_t len_threshold_;
            double absolute_barcode_threshold_;
            ScaffoldingUniqueEdgeStorage unique_storage_;

        public:
            TrivialTSLRExtensionChooser(const conj_graph_pack& gp, size_t len_threshold,
                double absolute_barcode_threshold, const ScaffoldingUniqueEdgeStorage& unique_storage) :
                    ExtensionChooser(gp.g),
                    bmapper_(gp.barcode_mapper),
                    len_threshold_(len_threshold),
                    absolute_barcode_threshold_(absolute_barcode_threshold),
                    unique_storage_(unique_storage) {
            }

            EdgeContainer Filter(const BidirectionalPath &path, const EdgeContainer &edges) const override {
                auto result = EdgeContainer();
                if (edges.size() == 0) {
                    return result;
                }
                //We might get single short edge as an input
                if (edges.size() == 1 && g_.length(edges.back().e_) < len_threshold_) {
                    result.push_back(edges.back());
                    return result;
                }
                //Find long unique edge earlier in path
                //FIXME code duplication
                bool long_single_edge_exists = false;
                EdgeId decisive_edge;
                for (int i = static_cast<int> (path.Size()) - 1; !long_single_edge_exists && i >= 0; --i) {
                    EdgeId current_edge = path[i];
                    if (unique_storage_.IsUnique(current_edge) &&
                            g_.length(current_edge) > len_threshold_) {
                        long_single_edge_exists = true;
                        decisive_edge = current_edge;
                    }
                }

                //Exclude decisive edge from the candidates
                auto edges_copy = edges;
                EraseEdge(edges_copy, decisive_edge);

                if (!long_single_edge_exists || edges_copy.size() == 0) {
                    if (edges_copy.size() == 0) {
                        DEBUG("Only decisive edge was found");
                    }
                    return result;
                }
                
                //Find edges with barcode score greater than some threshold
                std::vector <EdgeWithDistance> best_candidates;
                std::copy_if(edges_copy.begin(), edges_copy.end(), std::back_inserter(best_candidates),
                             [this, &decisive_edge](const EdgeWithDistance& edge) {
                                 return this->bmapper_->IntersectionSizeNormalizedBySecond(decisive_edge, edge.e_) >
                                        absolute_barcode_threshold_ && unique_storage_.IsUnique(edge.e_);
                             });
                if (best_candidates.size() == 1) {
                    result.push_back(best_candidates[0]);
                    return result;
                }
                if (best_candidates.size() == 0) {
                    return result;
                }


                //Check the difference between two best scores
                DEBUG("Several candidates found. Further filtering.");
                auto fittest_edge = *(std::max_element(edges_copy.begin(), edges_copy.end(),
                                                     [this, & decisive_edge](const EdgeWithDistance& edge1, const EdgeWithDistance& edge2) {
                                                         return this->bmapper_->IntersectionSizeNormalizedBySecond(decisive_edge, edge1.e_) <
                                                                this->bmapper_->IntersectionSizeNormalizedBySecond(decisive_edge, edge2.e_);
                                                     }));
                double best_score = bmapper_->IntersectionSizeNormalizedBySecond(decisive_edge, fittest_edge.e_);
                DEBUG("fittest edge " << fittest_edge.e_.int_id());
                DEBUG("score " << best_score);
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
                VERIFY(best_score >= second_best_score)
//                if (best_score - second_best_score < 0.03) { //todo configs
//                    DEBUG("Scores are too close, failed to select the best candidate.");
//                    return result;
//                }

                //Try to find topologically closest edge to resolve loops
                //FIXME This have nothing to do with barcodes. Need to be moved somewhere else.
                auto closest_edges = FindClosestEdge(best_candidates);
                if (closest_edges.size() != 1) {
                    DEBUG("Unable to find single topologically minimal edge.");
                    for (auto edge : best_candidates) {
                        result.push_back(edge);
                    }
                }
                else {
                    DEBUG("Found topologically minimal edge");
                    result.push_back(closest_edges.back());
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

            //make it more effective if needed
            vector<EdgeWithDistance> FindClosestEdge(const vector<EdgeWithDistance>& edges) const {
                vector <EdgeWithDistance> closest_edges;
                auto edges_iter = edges.begin();
                size_t path_len_bound = cfg::get().ts_res.topsort_bound;
                do {
                    auto edge = *edges_iter;
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
                        closest_edges.push_back(*edges_iter);
                    }
                    ++edges_iter;
                } while (edges_iter != edges.end());
                return closest_edges;
            }
            DECL_LOGGER("TslrExtensionChooser")
        };
} //tslr_resolver