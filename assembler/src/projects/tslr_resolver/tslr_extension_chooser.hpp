#pragma once

#include <common/modules/path_extend/extension_chooser.hpp>
#include <common/modules/path_extend/path_filter.hpp>
#include <common/modules/path_extend/path_extender.hpp>
#include <common/modules/path_extend/pe_resolver.hpp>
#include <common/modules/path_extend/path_extend_launch.hpp>
#include <assembly_graph/graph_support/scaff_supplementary.hpp>
#include <cmath>
#include "barcode_mapper.hpp"

using namespace path_extend;

namespace tslr_resolver {

        class ReadCloudExtensionChooser : public ExtensionChooser {

            shared_ptr<BarcodeMapper> bmapper_;
            size_t len_threshold_;
            double absolute_barcode_threshold_;
            size_t fragment_len_;
            ScaffoldingUniqueEdgeStorage unique_storage_;

        public:
            ReadCloudExtensionChooser(const conj_graph_pack& gp,
                                        size_t len_threshold,
                                        double absolute_barcode_threshold,
                                        size_t fragment_len,
                                        const ScaffoldingUniqueEdgeStorage& unique_storage) :
                    ExtensionChooser(gp.g),
                    bmapper_(gp.barcode_mapper),
                    len_threshold_(len_threshold),
                    absolute_barcode_threshold_(absolute_barcode_threshold),
                    fragment_len_(fragment_len),
                    unique_storage_(unique_storage) {
            }

            EdgeContainer Filter(const BidirectionalPath &path, const EdgeContainer &edges) const override {
                auto result = EdgeContainer();
                if (edges.size() == 0) {
                    return result;
                }
                //We might get a single short edge as an input
                if (edges.size() == 1 && g_.length(edges.back().e_) < len_threshold_) {
                    result.push_back(edges.back());
                    return result;
                }
                //Find last unique edge earlier in the path
                pair<EdgeId, int> last_unique = FindLastUniqueInPath(path, unique_storage_);
                bool long_single_edge_exists = false;
                EdgeId decisive_edge;
                if (last_unique.second != -1) {
                    long_single_edge_exists = true;
                    decisive_edge = last_unique.first;
                }

                //Exclude this edge from the candidates
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
                                 return this->bmapper_->GetIntersectionSizeNormalizedBySecond(decisive_edge, edge.e_) >
                                        absolute_barcode_threshold_ * GetGapCoefficient(edge.d_);
                             });
                if (best_candidates.size() == 1) {
                    result.push_back(best_candidates[0]);
                    DEBUG("Only one candidate passed threshold")
                    return result;
                }
                if (best_candidates.size() == 0) {
                    DEBUG("No candidates found")
                    return result;
                }
                if (bmapper_->GetTailBarcodeNumber(decisive_edge) < 4) {
                    DEBUG("Not enough barcodes on decisive edge")
                }


                //Check the difference between two best scores (currently inactive)
                DEBUG("Several candidates found. Further filtering.");
                auto best_edge = *(std::max_element(edges_copy.begin(), edges_copy.end(),
                                                     [this, & decisive_edge](const EdgeWithDistance& edge1,
                                                                             const EdgeWithDistance& edge2) {
                                                         return this->bmapper_->GetIntersectionSizeNormalizedBySecond(
                                                                 decisive_edge, edge1.e_) <
                                                                 this->bmapper_->GetIntersectionSizeNormalizedBySecond(
                                                                         decisive_edge, edge2.e_);
                                                     }));
                double best_score = bmapper_->GetIntersectionSizeNormalizedBySecond(decisive_edge, best_edge.e_);
                DEBUG("fittest edge " << best_edge.e_.int_id());
                DEBUG("score " << best_score);
                std::nth_element(edges_copy.begin(), edges_copy.begin() + 1, edges_copy.end(),
                                 [this, & decisive_edge](const EdgeWithDistance& edge1, const EdgeWithDistance& edge2) {
                                     return this->bmapper_->GetIntersectionSizeNormalizedBySecond(decisive_edge,
                                                                                                  edge1.e_) >
                                             this->bmapper_->GetIntersectionSizeNormalizedBySecond(decisive_edge,
                                                                                                   edge2.e_);
                                 });
                auto second_best_edge = edges_copy.at(1);
                double second_best_score = bmapper_->GetIntersectionSizeNormalizedBySecond(decisive_edge,
                                                                                           second_best_edge.e_);
                DEBUG("Second best edge " << second_best_edge.e_.int_id());
                DEBUG("second best score " << second_best_score);
                DEBUG(best_candidates.size() << " best candidates");
                VERIFY(best_score >= second_best_score)
//                if (best_score - second_best_score < 0.03) {
//                    DEBUG("Scores are too close, failed to select the best candidate.");
//                    return result;
//                }

                //Try to find topologically closest edge to resolve loops
                //fixme This have nothing to do with barcodes. Need to be moved elsewhere.
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

            //barcode threshold depends on gap between edges
            //Public version of this method
            static double GetGapCoefficient(int gap, size_t fragment_len) {
                VERIFY(gap <= (int)fragment_len)
                return static_cast<double>(fragment_len - gap) /
                       static_cast<double>(fragment_len);
            }

            //todo Does it really belong here?
            static std::pair<EdgeId, int> FindLastUniqueInPath(const BidirectionalPath& path,
                                                               const ScaffoldingUniqueEdgeStorage& storage) {
                for (int i =  (int)path.Size() - 1; i >= 0; --i) {
                    if (storage.IsUnique(path.At(i))) {
                        return std::make_pair(path.At(i), i);
                    }
                }
                return std::make_pair(EdgeId(0), -1);
            }

            static std::pair<EdgeId, int> FindFirstUniqueInPath(const BidirectionalPath& path,
                                                               const ScaffoldingUniqueEdgeStorage& storage) {
                for (int i = 0; i < (int)path.Size(); ++i) {
                    if (storage.IsUnique(path.At(i))) {
                        return std::make_pair(path.At(i), i);
                    }
                }
                return std::make_pair(EdgeId(0), -1);
            }

        private:

            //barcode threshold depends on gap between edges
            double GetGapCoefficient(int gap) const {
                VERIFY(gap <= (int)fragment_len_)
                return static_cast<double>(fragment_len_ - gap) /
                        static_cast<double>(fragment_len_);
            }

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