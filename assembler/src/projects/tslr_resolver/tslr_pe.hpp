#pragma once

#include <algorithms/path_extend/extension_chooser.hpp>
#include <algorithms/path_extend/path_filter.hpp>
#include <algorithms/path_extend/overlap_analysis.hpp>
#include <assembly_graph/graph_support/scaff_supplementary.hpp>
#include <cmath>
#include "barcode_mapper.hpp"
#include "tslr_visualizer.hpp"
#include "bounded_dijkstra.hpp"

#include <modules/algorithms/path_extend/path_extender.hpp>
#include <modules/algorithms/path_extend/pe_resolver.hpp>
#include <modules/algorithms/path_extend/path_extend_launch.hpp>

using namespace path_extend;

namespace tslr_resolver {

    struct res_statistics {
        size_t unaligned = 0;
        size_t true_pos = 0;
        size_t false_pos = 0;
        size_t negatives = 0;

        std::unordered_set <int64_t> true_positions;
        std::unordered_set <int64_t> false_positions;
        std::unordered_set <int64_t> neg_positions;

        void Serialize(const std::string filename) {
            std::ofstream fout;
            fout.open(filename);
            fout << "Unaligned: " << unaligned << endl;
            fout << "True positives: " << false_pos << endl <<
            "False positives: " << false_pos << endl <<
            "Negatives: " << negatives << endl;
            fout << "True positive positions" << endl;
            for (auto pos : false_positions) {
                fout << pos << ' ';
            }
            fout << endl;
            fout << "False positive positions" << endl;
            for (auto pos : false_positions) {
                fout << pos << ' ';
            }
            fout << endl;
            fout << "Negative positions" << endl;
            for (auto pos : neg_positions) {
                fout << pos << ' ';
            }
            fout << endl;
        }
    };

    class TrivialTSLRExtensionChooser : public ExtensionChooser {

        shared_ptr<BarcodeMapper> bmapper_;
        const EdgesPositionHandler<Graph>& edge_pos_;
        mutable res_statistics stats_;
        int reference_cov_;
        size_t len_threshold_;
        double relative_diff_threshold_;

    public:
        TrivialTSLRExtensionChooser(const conj_graph_pack& gp, const int& reference_cov, const size_t& len_threshold,
            const double& relative_diff_threshold) :
                ExtensionChooser(gp.g), bmapper_(gp.barcode_mapper), 
                edge_pos_(gp.edge_pos), stats_(), reference_cov_(reference_cov), 
                len_threshold_(len_threshold), relative_diff_threshold_(relative_diff_threshold) {
        }

        EdgeContainer Filter(const BidirectionalPath &path, const EdgeContainer &edges) const override {
            auto result = EdgeContainer();
            if (edges.size() == 0) {
                return result;
            }
            //Find long unique edge earlier in path
            int reference_cov = 20; 
            size_t len_threshold = 1000;
            double relative_diff_threshold = 0.05;
            bool long_single_edge_exists = false;
            EdgeId decisive_edge;
            for (int i = static_cast<int> (path.Size()) - 1; !long_single_edge_exists && i >= 0; --i) {
                EdgeId current_edge = path[i];
                if (g_.coverage(current_edge) < reference_cov * 1.1 &&
                    g_.length(current_edge) > len_threshold) {
                    long_single_edge_exists = true;
                    decisive_edge = current_edge;
                }
            }
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
            

            auto fittest_edge = *(std::max_element(edges_copy.begin(), edges_copy.end(),
                                                 [this, & decisive_edge](const EdgeWithDistance& edge1, const EdgeWithDistance& edge2) {
                                                     return this->bmapper_->IntersectionSizeRelative(decisive_edge, edge1.e_) <
                                                            this->bmapper_->IntersectionSizeRelative(decisive_edge, edge2.e_);
                                                 }));
            double best_score = bmapper_->IntersectionSizeRelative(decisive_edge, fittest_edge.e_);

            std::nth_element(edges_copy.begin(), edges_copy.begin() + 1, edges_copy.end(), 
                                                 [this, & decisive_edge](const EdgeWithDistance& edge1, const EdgeWithDistance& edge2) {
                                                     return this->bmapper_->IntersectionSizeRelative(decisive_edge, edge1.e_) >
                                                            this->bmapper_->IntersectionSizeRelative(decisive_edge, edge2.e_);
                                                 });
            auto second_best_edge = edges_copy.at(1);
            double second_best_score = bmapper_->IntersectionSizeRelative(decisive_edge, second_best_edge.e_);
            DEBUG("At edge " << path.Back().int_id());
            DEBUG("decisive edge " << decisive_edge.int_id());
            DEBUG("fittest edge " << fittest_edge.e_.int_id());
            DEBUG("score " << best_score);
            DEBUG("second best edge " << second_best_edge.e_.int_id());
            DEBUG("score " << second_best_score << endl);
            CheckAnswer(edge_pos_, decisive_edge, fittest_edge.e_, edges);
            result.push_back(fittest_edge);
            if (best_score - second_best_score < relative_diff_threshold)
                result.push_back(second_best_edge);
            return result;
        }

        void SerializeStats (const string& path) {
            stats_.Serialize(path);
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

        void CheckAnswer (const EdgesPositionHandler <Graph>& edge_pos, const EdgeId& decisive,
                          const EdgeId& candidate, const EdgeContainer& edges) const {
            if (edge_pos.GetEdgePositions(decisive).size() == 0) {
                stats_.unaligned++;
                return;
            }
            auto end_pos_contig = GetEndPos(edge_pos, decisive);
            string contig_name = end_pos_contig.first;
            int64_t initial_pos = end_pos_contig.second;
            EdgeId argmin = edges.begin()->e_;
            int64_t min = std::numeric_limits<int64_t>::max();
            bool is_forward_edge_found = false;
            for (auto it = edges.begin(); it != edges.end(); ++it) {
                auto edge = *it;
                if (edge_pos.GetEdgePositions(edge.e_).size() != 0) {
                    auto current_pos_contig = GetStartPos(edge_pos, edge.e_);
                    string current_contig_name = current_pos_contig.first;
                    int64_t current_pos = current_pos_contig.second;
                    if (contig_name == current_contig_name && current_pos > initial_pos) {
                        is_forward_edge_found = true;
                        if (initial_pos + min > current_pos) {
                            min = current_pos - initial_pos;
                            argmin = edge.e_;
                        }
                    }
                }
            }
            if (!is_forward_edge_found) {
                stats_.negatives++;
                stats_.neg_positions.insert(initial_pos);
                return;
            }
            if (argmin != candidate) {
                stats_.false_pos++;
                stats_.false_positions.insert(initial_pos);
            }
            else {
                stats_.true_pos++;
                stats_.true_positions.insert(initial_pos);
            }
        }

        std::pair <string, int64_t> GetEndPos (const EdgesPositionHandler <Graph>& edge_pos, const EdgeId& edge) const {
            EdgePosition endpos = edge_pos.GetEdgePositions(edge).back();
            return std::make_pair (endpos.contigId, endpos.mr.initial_range.end_pos);
        }

        std::pair <string, int64_t> GetStartPos (const EdgesPositionHandler <Graph>& edge_pos, const EdgeId& edge) const {
            EdgePosition startpos = edge_pos.GetEdgePositions(edge).front();
            return std::make_pair (startpos.contigId, startpos.mr.initial_range.end_pos);
        }

        DECL_LOGGER("TslrExtensionChooser")
    };

    class InconsistentTSLRExtender : public LoopDetectingPathExtender { //Traverse forward to find long edges

    protected:

        shared_ptr<ExtensionChooser> extensionChooser_;
        size_t distance_bound_;
        size_t edge_threshold_;

        void FindFollowingEdges(BidirectionalPath &path, ExtensionChooser::EdgeContainer *result) {
            result->clear();
            if (g_.OutgoingEdgeCount(g_.EdgeEnd(path.Back())) == 1) {
                result->push_back(EdgeWithDistance(*(g_.OutgoingEdges(g_.EdgeEnd(path.Back())).begin()), 0));
                return;
            }
            auto dij = LengthDijkstra<Graph>::CreateLengthBoundedDijkstra(g_, distance_bound_, edge_threshold_);
            dij.Run(g_.EdgeEnd(path.Back()));
            auto processed_vertices = dij.ProcessedVertices();
            std::set<EdgeId> long_edges;
            for (auto vertex : processed_vertices) {
                for (auto edge : g_.OutgoingEdges(vertex)) {
                    if (g_.length(edge) >= edge_threshold_) {
                        long_edges.insert(edge);
                    }
                }
            }
            result->reserve(long_edges.size());
            for (auto edge : long_edges) {
                result->push_back(EdgeWithDistance(edge, dij.GetDistance(g_.EdgeStart(edge))));
            }
        }


    public:

        InconsistentTSLRExtender(const conj_graph_pack &gp, const GraphCoverageMap &cov_map, shared_ptr<ExtensionChooser> ec,
                                       size_t is, size_t max_loops, bool investigate_short_loops, bool use_short_loop_cov_resolver, 
                                       size_t distance_bound, size_t edge_threshold)
                :
                LoopDetectingPathExtender(gp, cov_map, max_loops, investigate_short_loops, use_short_loop_cov_resolver,
                                          is),
                extensionChooser_(ec), distance_bound_(distance_bound), edge_threshold_(edge_threshold) {
        }

        std::shared_ptr<ExtensionChooser> GetExtensionChooser() const {
            return extensionChooser_;
        }

        bool CanInvestigateShortLoop() const override {
            return extensionChooser_->WeightCounterBased();
        }

        bool ResolveShortLoopByCov(BidirectionalPath &path) override {
            return path.Size() < 1;
        }

        bool ResolveShortLoopByPI(BidirectionalPath &path) override {
            return path.Size() < 1;
        }

        bool MakeSimpleGrowStep(BidirectionalPath &path, PathContainer *paths_storage) override {
            ExtensionChooser::EdgeContainer candidates;
            return FilterCandidates(path, candidates) and AddCandidates(path, paths_storage, candidates);
        }

    protected:
        virtual bool FilterCandidates(BidirectionalPath &path, ExtensionChooser::EdgeContainer &candidates) {
            if (path.Size() == 0) {
                return false;
            }
            DEBUG("Simple grow step");
            path.Print();
            DEBUG("Starting at vertex " << g_.EdgeEnd(path.Back()));
            FindFollowingEdges(path, &candidates);
            DEBUG("found candidates");
            DEBUG(candidates.size())
            if (candidates.size() == 1) {
                LoopDetector loop_detector(&path, cov_map_);
                if (!investigate_short_loops_ &&
                    (loop_detector.EdgeInShortLoop(path.Back()) or loop_detector.EdgeInShortLoop(candidates.back().e_))
                    && extensionChooser_->WeightCounterBased()) {
                    return false;
                }
            }
            DEBUG("more filtering");
            candidates = extensionChooser_->Filter(path, candidates);
            DEBUG("filtered candidates");
            DEBUG(candidates.size())
            return true;
        }

        virtual bool AddCandidates(BidirectionalPath &path, PathContainer * /*paths_storage*/,
                                   ExtensionChooser::EdgeContainer &candidates) {
            if (candidates.size() != 1) {
                if (candidates.size() > 1) {
                    DEBUG("Too many candidates, false");
                }
                if (candidates.size() == 0) {
                    DEBUG("No candidates found, false");
                }
                DEBUG("Final(?) path length: " << path.Length());
                return false;
            }

            LoopDetector loop_detector(&path, cov_map_);
            DEBUG("loop detector");
            if (!investigate_short_loops_ &&
                (loop_detector.EdgeInShortLoop(path.Back()) or loop_detector.EdgeInShortLoop(candidates.back().e_))
                && extensionChooser_->WeightCounterBased()) {
                return false;
            }
            DEBUG("push");
            EdgeId eid = candidates.back().e_;
            //In 2015 modes when trying to use already used unique edge, it is not added and path growing stops.
            //That allows us to avoid overlap removal hacks used earlier.
            if (used_storage_->UniqueCheckEnabled()) {
                if (used_storage_->IsUsedAndUnique(eid)) {
                    return false;
                } else {
                    used_storage_->insert(eid);
                }
            }
            path.PushBack(eid, candidates.back().d_);
            DEBUG("Push done, true");
            DEBUG("Path length: " << path.Length());
            return true;
        }

    protected:
        DECL_LOGGER("InconsistentExtender")

    };

    void LaunchBarcodePE (conj_graph_pack &gp) {
        path_extend::PathExtendParamsContainer params(cfg::get().pe_params,
                                                      cfg::get().output_dir,
                                                      "final_contigs_tslr",
                                                      "scaffolds_tslr",
                                                      cfg::get().mode,
                                                      cfg::get().uneven_depth,
                                                      cfg::get().avoid_rc_connections,
                                                      cfg::get().use_scaffolder);

        DefaultContigCorrector<ConjugateDeBruijnGraph> corrector(gp.g);
        DefaultContigConstructor<ConjugateDeBruijnGraph> constructor(gp.g, corrector);
        ContigWriter writer(gp.g, constructor, gp.components, params.mode == config::pipeline_type::plasmid);
        GraphCoverageMap clone_map(gp.g);
        const pe_config::ParamSetT &pset = params.pset;
        bool detect_repeats_online = false;
        int reference_cov = 20;
        size_t len_threshold = 1000;
        double relative_diff_threshold = 0.05;

        auto extension = make_shared<TrivialTSLRExtensionChooser>(gp, reference_cov, len_threshold, relative_diff_threshold);
        size_t distance_bound = 10000; //TODO configs
        size_t edge_threshold = 1000;
        auto tslr_extender = make_shared<InconsistentTSLRExtender>(gp, clone_map,
                                                             extension,
                                                             2500 /*insert size*/,
                                                             0 /*max loops*/,
                                                             false, /*investigate short loops*/
                                                             false /*use short loop coverage resolver*/,
                                                             distance_bound,
                                                             edge_threshold);
        vector <shared_ptr<PathExtender> > all_libs;
        all_libs.push_back(tslr_extender);
        size_t max_is_right_quantile = gp.g.k() + 10000;
        ScaffoldingUniqueEdgeStorage main_unique_storage;
        shared_ptr<CompositeExtender> mainPE = make_shared<CompositeExtender>
                (gp.g, clone_map, all_libs,
                 main_unique_storage,
                 max_is_right_quantile,
                 pset.extension_options.max_repeat_length,
                 detect_repeats_online);

        PathExtendResolver resolver(gp.g);
        auto seeds = resolver.makeSimpleSeeds();
        seeds.SortByLength();
        auto last_paths = resolver.extendSeeds(seeds, *mainPE);

        extension->SerializeStats(cfg::get().output_dir + "stats");
        size_t min_edge_len = 100;
        FinalizePaths(params, last_paths, clone_map, min_edge_len, max_is_right_quantile);

        debruijn_graph::GenomeConsistenceChecker genome_checker (gp, main_unique_storage, 1000, 0.2);
        DebugOutputPaths(gp, params, last_paths, "final_tslr_paths");

        writer.OutputPaths(last_paths, params.output_dir + params.contigs_name);
        if (gp.genome.size() > 0)
            CountMisassembliesWithReference(genome_checker, last_paths);
    };

} //tslr_resolver
