#pragma once

#include "tslr_extension_chooser.hpp"
#include "tslr_visualizer.hpp"
#include "bounded_dijkstra.hpp"
#include "algorithms/dijkstra/dijkstra_helper.hpp"

using namespace path_extend;

namespace tslr_resolver {

    class InconsistentTSLRExtender : public LoopDetectingPathExtender { //Traverse forward to find long edges

    protected:

        shared_ptr<ExtensionChooser> extensionChooser_;
        size_t distance_bound_;
        double barcode_threshold_;
        size_t edge_threshold_;
        size_t reference_cov_;
        shared_ptr<tslr_resolver::BarcodeMapper> mapper_;


        void FindFollowingEdges(BidirectionalPath &path, ExtensionChooser::EdgeContainer *result) {
            result->clear();
            if (g_.OutgoingEdgeCount(g_.EdgeEnd(path.Back())) == 1) {
                result->push_back(EdgeWithDistance(*(g_.OutgoingEdges(g_.EdgeEnd(path.Back())).begin()), 0));
                return;
            }

            //find long unique edge earlier in path
            bool long_single_edge_exists = false;
            EdgeId decisive_edge;
            for (int i = static_cast<int> (path.Size()) - 1; !long_single_edge_exists && i >= 0; --i) {
                EdgeId current_edge = path[i];
                if (IsEdgeUnique(g_, current_edge) &&
                        g_.length(current_edge) > edge_threshold_) {
                    long_single_edge_exists = true;
                    decisive_edge = current_edge;
                }
            }

            DEBUG("At edge " << path.Back().int_id());
            DEBUG("decisive edge " << decisive_edge.int_id());

            if (!long_single_edge_exists) {
                DEBUG("Couldn't find single long edge");
                return;
            }

            DEBUG("decisive edge barcodes: " <<  mapper_->GetSizeTails(decisive_edge));

            //find reliable unique edges further in graph
            vector <EdgeId> candidates;
            auto put_checker = BarcodePutChecker<Graph>(g_, edge_threshold_, barcode_threshold_, 
                mapper_, decisive_edge, candidates);
            auto dij = BarcodeDijkstra<Graph>::CreateBarcodeBoundedDijkstra(g_, distance_bound_, put_checker);
            dij.Run(g_.EdgeEnd(path.Back()));
            result->reserve(candidates.size());
            for (auto edge : candidates) {
                result->push_back(EdgeWithDistance(edge, dij.GetDistance(g_.EdgeStart(edge))));
            }
        }


    public:

        InconsistentTSLRExtender(const conj_graph_pack &gp, const GraphCoverageMap &cov_map, shared_ptr<ExtensionChooser> ec,
                                       size_t is, size_t max_loops, bool investigate_short_loops, bool use_short_loop_cov_resolver, 
                                       size_t distance_bound, double barcode_threshold, size_t edge_threshold)
                :
                LoopDetectingPathExtender(gp, cov_map, max_loops, investigate_short_loops, use_short_loop_cov_resolver,
                                          is),
                extensionChooser_(ec), distance_bound_(distance_bound), barcode_threshold_(barcode_threshold), 
                    edge_threshold_(edge_threshold), mapper_(gp.barcode_mapper) {
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

    PathContainer FilterByLength(const PathContainer& seeds, size_t len_threshold) {
        PathContainer result;
        for (auto it = seeds.begin(); it != seeds.end(); ++it) {
            if (it->first->Length() > len_threshold) {
                result.AddPair(it->first, it->second);
            }
        }
        return result;
    }

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
        GraphCoverageMap cover_map(gp.g);
        const pe_config::ParamSetT &pset = params.pset;
        bool detect_repeats_online = false;

        PathExtendResolver resolver(gp.g);
        auto min_unique_length = pset.scaffolding2015.min_unique_length;
        auto unique_variaton = pset.scaffolding2015.unique_coverage_variation;
        ScaffoldingUniqueEdgeStorage main_unique_storage;
        auto dataset_info = cfg::get().ds;
        main_unique_storage = FillUniqueEdgeStorage(gp, dataset_info,
                                                    min_unique_length,
                                                    unique_variaton,
                                                    false);

        //mp extender
        INFO("SUBSTAGE = paired-end libraries")
        PathExtendStage exspander_stage = PathExtendStage::PEStage;
        vector<shared_ptr<PathExtender> > all_libs =
            MakeAllExtenders(exspander_stage, dataset_info, params, gp, cover_map, main_unique_storage);
        size_t max_is_right_quantile = max(FindOverlapLenForStage(exspander_stage, cfg::get().ds), gp.g.k() + 100);
        size_t min_edge_len = 100;

        shared_ptr<CompositeExtender> mainPE = make_shared<CompositeExtender>(gp.g, cover_map, all_libs,
                                                                              main_unique_storage,
                                                                              max_is_right_quantile,
                                                                              pset.extension_options.max_repeat_length,
                                                                              detect_repeats_online);
        auto seeds = resolver.makeSimpleSeeds();
        seeds.SortByLength();
        auto paths = resolver.extendSeeds(seeds, *mainPE);
        paths.SortByLength();
        FinalizePaths(params, paths, cover_map, min_edge_len, max_is_right_quantile);
        seeds.DeleteAllPaths();
        all_libs.clear();

        //tslr extender
        INFO("SUBSTAGE = TSLR Resolver")
        auto tslr_resolver_params = cfg::get().ts_res; //TODO: move params reading somewhere
        int reference_cov = tslr_resolver_params.reference_cov;
        size_t len_threshold = tslr_resolver_params.len_threshold;
        double relative_diff_threshold = tslr_resolver_params.diff_threshold;
        size_t distance_bound = tslr_resolver_params.distance_bound;
        double abs_threshold = tslr_resolver_params.abs_threshold;

        max_is_right_quantile = gp.g.k() + 10000;
        auto extension = make_shared<TrivialTSLRExtensionChooser>(gp, reference_cov, len_threshold, relative_diff_threshold);
        auto tslr_extender = make_shared<InconsistentTSLRExtender>(gp, cover_map,
                                                             extension,
                                                             2500 /*insert size*/,
                                                             0 /*max loops*/,
                                                             false, /*investigate short loops*/
                                                             false /*use short loop coverage resolver*/,
                                                             distance_bound,
                                                             abs_threshold,
                                                             len_threshold);
        all_libs.push_back(tslr_extender);
        shared_ptr<CompositeExtender> tslrPE = make_shared<CompositeExtender>
                (gp.g, cover_map, all_libs,
                 main_unique_storage,
                 max_is_right_quantile,
                 pset.extension_options.max_repeat_length,
                 detect_repeats_online);

        auto filtered_paths = FilterByLength(paths, len_threshold);
        auto tslr_paths = resolver.extendSeeds(filtered_paths, *tslrPE);

        FinalizePaths(params, tslr_paths, cover_map, min_edge_len, max_is_right_quantile);

        debruijn_graph::GenomeConsistenceChecker genome_checker (gp, main_unique_storage, 1000, 0.2);
        DebugOutputPaths(gp, params, tslr_paths, "final_tslr_paths");

        writer.OutputPaths(tslr_paths, params.output_dir + params.contigs_name);
        if (gp.genome.size() > 0)
            CountMisassembliesWithReference(genome_checker, tslr_paths);
    };

} //tslr_resolver
