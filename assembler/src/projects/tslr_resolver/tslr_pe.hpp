#pragma once

#include <algorithms/path_extend/extension_chooser.hpp>
#include <algorithms/path_extend/path_filter.hpp>
#include <algorithms/path_extend/overlap_analysis.hpp>
#include <assembly_graph/graph_support/scaff_supplementary.hpp>
#include <cmath>
#include "barcode_mapper.hpp"

#include <modules/algorithms/path_extend/path_extender.hpp>
#include <modules/algorithms/path_extend/pe_resolver.hpp>
#include <modules/algorithms/path_extend/path_extend_launch.hpp>

using namespace path_extend;

namespace tslr_resolver {
    class TrivialTSLRExtensionChooser : public ExtensionChooser {

        const BarcodeMapper& bmapper_;

    public:
        TrivialTSLRExtensionChooser(Graph &g, const BarcodeMapper& bmapper) :
                ExtensionChooser(g), bmapper_(bmapper) {
        }

        EdgeContainer Filter(const BidirectionalPath &path, const EdgeContainer &edges) const override {
            auto result = EdgeContainer();
            if (edges.size() == 1) {
            	DEBUG("on branch 1");
                return edges;
            }
            int reference_cov = 20;
            size_t len_threshold = 1000;
            bool long_single_edge_exists = false;
            EdgeId decisive_edge;
            for (int i = static_cast<int> (path.Size()) - 1; !long_single_edge_exists && i >= 0; --i) {
                EdgeId current_edge = path[i];
                if (g_.coverage(current_edge) < reference_cov * 1.5 &&
                    g_.length(current_edge) < len_threshold) {
                    long_single_edge_exists = true;
                    decisive_edge = current_edge;
                }
            }
            if (!long_single_edge_exists || edges.size() == 0) {
            	DEBUG("on branch 2");
                return result;
            }
            DEBUG("On branch 3");
            auto fittest_edge = std::max_element(edges.begin(), edges.end(),
                                                 [this, & decisive_edge](const EdgeWithDistance& edge1, const EdgeWithDistance& edge2) {
                                                     return this->bmapper_.IntersectionSize(edge1.e_, decisive_edge) <
                                                            this->bmapper_.IntersectionSize(edge2.e_, decisive_edge);
                                                 });
            result.push_back(*fittest_edge);
            return result;
        }
        DECL_LOGGER("TslrExtensionChooser")
    };

    class SimpleTSLRExtender : public LoopDetectingPathExtender { //Same as SimpleExtender, but with removed loop checks

    protected:

        shared_ptr<ExtensionChooser> extensionChooser_;

        void FindFollowingEdges(BidirectionalPath &path, ExtensionChooser::EdgeContainer *result) {
            DEBUG("Looking for the following edges")
            result->clear();
            vector<EdgeId> edges;
            DEBUG("Pushing back")
            push_back_all(edges, g_.OutgoingEdges(g_.EdgeEnd(path.Back())));
            result->reserve(edges.size());
            for (auto iter = edges.begin(); iter != edges.end(); ++iter) {
                DEBUG("Adding edge w distance " << g_.int_id(*iter));
                result->push_back(EdgeWithDistance(*iter, 0));
            }
            DEBUG("Following edges found");
        }


    public:

        SimpleTSLRExtender(const conj_graph_pack &gp, const GraphCoverageMap &cov_map, shared_ptr<ExtensionChooser> ec,
                           size_t is, size_t max_loops, bool investigate_short_loops, bool use_short_loop_cov_resolver)
                :
                LoopDetectingPathExtender(gp, cov_map, max_loops, investigate_short_loops, use_short_loop_cov_resolver,
                                          is),
                extensionChooser_(ec) {
        }

        std::shared_ptr<ExtensionChooser> GetExtensionChooser() const {
            return extensionChooser_;
        }

        bool CanInvestigateShortLoop() const override {
            return extensionChooser_->WeightCounterBased();
        }

        bool ResolveShortLoopByCov(BidirectionalPath &path) override {
            return path.Size() < 1; //to avoid warning, will be removed later
        }

        bool ResolveShortLoopByPI(BidirectionalPath &path) override {
            return path.Size() < 1; //to avoid warning, will be removed later
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
            if (candidates.size() != 1)
                return false;

            LoopDetector loop_detector(&path, cov_map_);
            DEBUG("loop detecor");
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
            DEBUG("push done");
            return true;
        }

    protected:
        DECL_LOGGER("SimpleExtender")

    };

    void LaunchBarcodePE (conj_graph_pack &gp) {
        path_extend::PathExtendParamsContainer params(cfg::get().pe_params,
                                                      cfg::get().output_dir,
                                                      "final_contigs",
                                                      "scaffolds",
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
        auto extension = make_shared<TrivialTSLRExtensionChooser>(gp.g, gp.barcode_mapper);
        auto tslr_extender = make_shared<SimpleTSLRExtender>(gp, clone_map,
                                                             extension,
                                                             0 /*insert size*/,
                                                             0 /*max loops*/,
                                                             false, /*investigate short loops*/
                                                             false /*use short loop coverage resolver*/);
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

//extend pe + long reads
        PathExtendResolver resolver(gp.g);
        auto seeds = resolver.makeSimpleSeeds();
        seeds.SortByLength();
        INFO("Growing paths using paired-end and long single reads");
        auto last_paths = resolver.extendSeeds(seeds, *mainPE);

        debruijn_graph::GenomeConsistenceChecker genome_checker (gp, main_unique_storage, 1000, 0.2);
        DebugOutputPaths(gp, params, last_paths, "final_tslr_paths");

        writer.OutputPaths(last_paths, params.output_dir + params.contigs_name);
        if (gp.genome.size() > 0)
            CountMisassembliesWithReference(genome_checker, last_paths);
    };

} //tslr_resolver
