//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "simplification.hpp"

#include "assembly_graph/core/basic_graph_stats.hpp"
#include "assembly_graph/graph_support/graph_processing_algorithm.hpp"
#include "assembly_graph/stats/picture_dump.hpp"

#include "modules/simplification/cleaner.hpp"

#include "stages/simplification_pipeline/simplification_settings.hpp"
#include "stages/simplification_pipeline/graph_simplification.hpp"
#include "stages/simplification_pipeline/single_cell_simplification.hpp"
#include "stages/simplification_pipeline/rna_simplification.hpp"

#include "pipeline/genomic_info.hpp"

#include "utils/perf/timetracer.hpp"

namespace debruijn_graph {

using namespace debruijn::simplification;
using namespace config;

SimplifInfoContainer CreateInfoContainer(const GraphPack &gp) {
    SimplifInfoContainer info_container(cfg::get().mode);
    info_container.set_read_length(cfg::get().ds.RL)
            .set_main_iteration(cfg::get().main_iteration)
            .set_chunk_cnt(5 * cfg::get().max_threads);

    //0 if model didn't converge
    //todo take max with trusted_bound
    const auto &ginfo = gp.get<GenomicInfo>();
    info_container.set_detected_coverage_bound(ginfo.ec_bound());
    if (!cfg::get().uneven_depth) {
        info_container.set_detected_mean_coverage(ginfo.estimated_mean());
    }

    return info_container;
}

class GraphSimplifier {
    typedef std::function<void(EdgeId)> HandlerF;
    typedef SmartEdgeSet<std::unordered_set<EdgeId>, Graph> RestrictedEdgeSet;

    GraphPack& gp_;
    Graph& g_;
    SimplifInfoContainer info_container_;
    const debruijn_config::simplification simplif_cfg_;
    HandlerF removal_handler_;
    stats::detail_info_printer& printer_;
    CountingCallback<Graph> cnt_callback_;

    bool PerformInitCleaning() {
        if (simplif_cfg_.init_clean.early_it_only && info_container_.main_iteration()) {
            INFO("Most init cleaning disabled on main iteration");
            return false;
        }
        if (math::ge(simplif_cfg_.init_clean.activation_cov, 0.)
                && math::ls(info_container_.detected_mean_coverage(), simplif_cfg_.init_clean.activation_cov)) {
            INFO("Most init cleaning disabled since detected mean " << info_container_.detected_mean_coverage()
                 << " was less than activation coverage " << simplif_cfg_.init_clean.activation_cov);
            return false;
        }

        return true;
    }

    bool FinalRemoveErroneousEdges() {
        CountingCallback<Graph> counting_callback;
        HandlerF handler = [&] (EdgeId e) {
            if (removal_handler_)
                removal_handler_(e);
            counting_callback.HandleDelete(e);
        };

        INFO("Topology-based removal procedures");
        bool changed = TopologyRemoveErroneousEdges(g_, simplif_cfg_.tec, handler);
        changed |= TopologyReliabilityRemoveErroneousEdges(g_, simplif_cfg_.trec, handler);
        changed |= RemoveThorns(g_, simplif_cfg_.isec, removal_handler_);
        changed |= MultiplicityCountingRemoveErroneousEdges(g_, simplif_cfg_.tec, handler);
        changed |= MaxFlowRemoveErroneousEdges(g_, simplif_cfg_.mfec, handler);
        counting_callback.Report();
        return changed;
    }

public:
    GraphSimplifier(GraphPack &gp, const SimplifInfoContainer& info_container,
                    const debruijn_config::simplification& simplif_cfg,
                    const std::function<void(EdgeId)>& removal_handler,
                    stats::detail_info_printer& printer)
            : gp_(gp),
              g_(gp.get_mutable<Graph>()),
              info_container_(info_container),
              simplif_cfg_(simplif_cfg),
              removal_handler_(removal_handler),
              printer_(printer) {
    }

    void InitialCleaningRNA() {
        INFO("PROCEDURE == Initial cleaning");
        TIME_TRACE_SCOPE("Initial cleaing (RNA)");

        printer_(info_printer_pos::before_raw_simplification);

        CompositeAlgorithm<Graph> algo(g_);

        algo.AddAlgo(LowComplexityShortEdgeRemoverInstance(g_, removal_handler_, info_container_.chunk_cnt()),
                     "poly A/T short edge remover");
        algo.AddAlgo(LowComplexityTipClipperInstance(g_, removal_handler_, info_container_.chunk_cnt()),
                     "poly A/T tip clipper");

        //Currently need force_primary_launch = true for correct behavior of looped TipClipper
        algo.Run(/*force primary launch*/true);

        const auto &flanking_cov = gp_.get<FlankingCoverage<Graph>>();
        RemoveHiddenLoopEC(g_, flanking_cov, info_container_.detected_coverage_bound(),
                           simplif_cfg_.her.relative_threshold, removal_handler_);
    }

    void InitialCleaning() {
        INFO("PROCEDURE == Initial cleaning");
        TIME_TRACE_SCOPE("Initial cleaing");

        printer_(info_printer_pos::before_raw_simplification);

        CompositeAlgorithm<Graph> algo(g_);
        const auto &flanking_cov = gp_.get<FlankingCoverage<Graph>>();

        algo.AddAlgo(
                SelfConjugateEdgeRemoverInstance(g_,
                                                 simplif_cfg_.init_clean.self_conj_condition,
                                                 info_container_, removal_handler_),
                "Self conjugate edge remover");

        if (PerformInitCleaning()) {
            algo.AddAlgo(
                    TipClipperInstance(g_,
                                       debruijn_config::simplification::tip_clipper(simplif_cfg_.init_clean.tip_condition),
                                       info_container_,
                                       removal_handler_),
                    "Initial tip clipper");

            algo.AddAlgo(
                    ECRemoverInstance(g_,
                                      debruijn_config::simplification::erroneous_connections_remover(simplif_cfg_.init_clean.ec_condition),
                                      info_container_,
                                      removal_handler_),
                    "Initial ec remover");

            algo.AddAlgo(
                    LowFlankDisconnectorInstance(g_, flanking_cov,
                                                 simplif_cfg_.init_clean.disconnect_flank_cov, info_container_,
                                                 removal_handler_),
                    "Disconnecting edges with low flanking coverage");

            algo.AddAlgo(
                    IsolatedEdgeRemoverInstance(g_,
                                                simplif_cfg_.init_clean.ier,
                                                info_container_, removal_handler_),
                    "Initial isolated edge remover");
        }

        //Currently need force_primary_launch = true for correct behavior of looped TipClipper
        algo.Run(/*force primary launch*/true);
    }

    void PostSimplification() {
        using namespace omnigraph;
        using namespace func;
        INFO("PROCEDURE == Post simplification");
        TIME_TRACE_SCOPE("Post simplification");

        printer_(info_printer_pos::before_post_simplification);

        //auto colorer = debruijn_graph::DefaultGPColorer(gp_);
        //visualization::graph_labeler::DefaultLabeler<Graph> labeler(g_, gp_.edge_pos);

        //    gp.ClearQuality();
        //    gp.FillQuality();
        //    QualityEdgeLocalityPrintingRH<Graph> qual_removal_handler(gp.g, gp.edge_qual, labeler, colorer,
        //                                   cfg::get().output_dir + "pictures/colored_edges_deleted/");
        //
        //    //positive quality edges removed (folder colored_edges_deleted)
        //    std::function<void(EdgeId)> qual_removal_handler_f = boost::bind(
        //            //            &QualityLoggingRemovalHandler<Graph>::HandleDelete,
        //            &QualityEdgeLocalityPrintingRH<Graph>::HandleDelete,
        //            boost::ref(qual_removal_handler), _1);

        //visualization::visualization_utils::LocalityPrintingRH<Graph> drawing_handler(gp_.g, labeler, colorer, "/home/snurk/pics");
        //auto printing_handler=[&] (EdgeId e) {
        //    std::cout << "Edge:" << g_.str(e) << "; cov: " << g_.coverage(e) << "; start " << g_.str(g_.EdgeStart(e)) << "; end " << g_.str(g_.EdgeEnd(e)) << std::endl;
        //};
        //auto extensive_handler = [&] (EdgeId e) {removal_handler_(e) ; printing_handler(e); drawing_handler.HandleDelete(e);};
        bool use_restricted = gp_.count<SmartEdgeSet<std::unordered_set<EdgeId>, Graph>>();
        if (use_restricted) {
            DEBUG("RestrictedEdgeSet is present");
        }

        typedef std::function<bool(EdgeId edge, const std::vector<EdgeId>& path)> BulgeCallbackF;
        BulgeCallbackF bulge_callback_f = nullptr;
        if (use_restricted) {
            bulge_callback_f = [this](EdgeId e, const std::vector<EdgeId>&) {
                                   DEBUG("Checking if " << e << " is in restricted edge set");
                                   return gp_.get<SmartEdgeSet<std::unordered_set<EdgeId>, Graph>>().count(e) > 0;
                               };
        }

        typename ComponentRemover<Graph>::HandlerF set_removal_handler_f;
        if (removal_handler_) {
            set_removal_handler_f = [=](const std::set<EdgeId> &edges) {
                std::for_each(edges.begin(), edges.end(), removal_handler_);
            };
        }

        CompositeAlgorithm<Graph> algo(g_);
        const auto &flanking_cov = gp_.get<FlankingCoverage<Graph>>();

        if (simplif_cfg_.topology_simplif_enabled && info_container_.main_iteration()) {
            algo.AddAlgo<AdapterAlgorithm<Graph>>("Final removal of erroneous edges", g_, [&]() { return FinalRemoveErroneousEdges(); });
            algo.AddAlgo(
                    TopologyTipClipperInstance(g_, simplif_cfg_.ttc,
                                               info_container_, removal_handler_),
                    "Topology tip clipper");
        }

        algo.AddAlgo(
                RelativeECRemoverInstance(g_, simplif_cfg_.rcec, info_container_, removal_handler_),
                "Relative coverage erroneous connection remover");

        algo.AddAlgo(
                RelativeCoverageComponentRemoverInstance(g_, flanking_cov,
                                                         simplif_cfg_.rcc, info_container_, set_removal_handler_f),
                "Relative coverage component remover");

        algo.AddAlgo(
                RelativelyLowCoverageDisconnectorInstance(g_,
                                                          simplif_cfg_.red, info_container_,
                                                          removal_handler_),
                "Disconnecting edges with relatively low coverage");

        algo.AddAlgo(
                ComplexTipClipperInstance(g_, simplif_cfg_.complex_tc, info_container_, set_removal_handler_f),
                "Complex tip clipper");

        algo.AddAlgo(
            ComplexBRInstance(g_, simplif_cfg_.cbr, use_restricted ? &gp_.get<SmartEdgeSet<std::unordered_set<EdgeId>, Graph>>() : nullptr, info_container_),
                "Complex bulge remover");

        algo.AddAlgo(
                TipClipperInstance(g_, simplif_cfg_.tc,
                                   info_container_, removal_handler_),
                "Tip clipper");

        algo.AddAlgo(
                TipClipperInstance(g_, simplif_cfg_.final_tc,
                                   info_container_, removal_handler_),
                "Final tip clipper");

        algo.AddAlgo(
                BRInstance(g_, simplif_cfg_.br,
                           info_container_, bulge_callback_f, removal_handler_),
                "Bulge remover");

        algo.AddAlgo(
                BRInstance(g_, simplif_cfg_.final_br,
                           info_container_, bulge_callback_f, removal_handler_),
                "Final bulge remover");
        algo.AddAlgo(
                BRInstance(g_, simplif_cfg_.subspecies_br,
                           info_container_, bulge_callback_f, removal_handler_),
                "Subspecies bulge remover");

        //TODO need better configuration
        if (config::PipelineHelper::IsMetagenomicPipeline(info_container_.mode())) {
            EdgePredicate<Graph> meta_thorn_condition
                    = And(LengthUpperBound<Graph>(g_, LengthThresholdFinder::MaxErroneousConnectionLength(
                                                                           g_.k(), simplif_cfg_.isec.max_ec_length_coefficient)),

                      And([&] (EdgeId e) {
                              //todo configure!
                              return simplification::relative_coverage::
                                         RelativeCoverageHelper<Graph>(g_, flanking_cov, 2).AnyHighlyCoveredOnFourSides(e);
                          },

                      And(UniqueIncomingPathLengthLowerBound(g_, simplif_cfg_.isec.uniqueness_length),
                          //todo configure!
                          TopologicalThornCondition<Graph>(g_, simplif_cfg_.isec.span_distance, /*max edge cnt*/5))));

            algo.AddAlgo<ParallelEdgeRemovingAlgorithm<Graph>>("Thorn remover (meta)", g_, meta_thorn_condition,
                                                               info_container_.chunk_cnt(),
                                                               removal_handler_);
        }

        // TODO: better configuration
        if (info_container_.mode() == config::pipeline_type::rna)
            algo.AddAlgo(LowComplexityTipClipperInstance(g_, removal_handler_, info_container_.chunk_cnt()), "AT Tips");

        algo.AddAlgo(
                LowCoverageEdgeRemoverInstance(g_,
                                               simplif_cfg_.lcer,
                                               info_container_),
                "Removing edges with low coverage");

        AlgorithmRunningHelper<Graph>::LoopedRunPrimaryOpening(algo,
                /*first primary iteration cnt*/2, /*max it count*/10);

        //TODO make part of cycle?
        RemoveHiddenEC(g_, flanking_cov, simplif_cfg_.her, info_container_, removal_handler_);

        //TODO better configuration
        if (config::PipelineHelper::IsMetagenomicPipeline(info_container_.mode())) {
            VERIFY(math::ls(simplif_cfg_.her.unreliability_threshold, 0.));
            MetaHiddenECRemover<Graph> algo(g_, info_container_.chunk_cnt(), flanking_cov,
                                            simplif_cfg_.her.uniqueness_length,
                                            simplif_cfg_.her.relative_threshold,
                                            removal_handler_);
            INFO("Running Hidden EC remover (meta)");
            AlgorithmRunningHelper<Graph>::LoopedRun(algo);
        }

        INFO("Disrupting self-conjugate edges");
        SelfConjugateDisruptor<Graph>(g_, cfg::get().max_repeat_length, removal_handler_).Run();

        AlgorithmRunningHelper<Graph>::RunAlgo(IsolatedEdgeRemoverInstance(g_, cfg::get().simp.ier,
                                                                           info_container_,
                                                                           removal_handler_),
                                               "Removing isolated edges");

        double low_threshold = gp_.get<GenomicInfo>().trusted_bound();
        if (math::gr(low_threshold, 0.0)) {
            INFO("Removing all the edges having coverage " << low_threshold << " and less");
            ParallelEdgeRemovingAlgorithm<Graph, CoverageComparator<Graph>>
                    cov_cleaner(g_,
                                CoverageUpperBound<Graph>(g_, low_threshold),
                                info_container_.chunk_cnt(),
                                removal_handler_,
                                /*canonical_only*/true,
                                CoverageComparator<Graph>(g_));
            cov_cleaner.Run();
        }

        printer_(info_printer_pos::final_simplified);
    }

    void SimplifyGraph() {
        TIME_TRACE_SCOPE("Graph simplification");

        bool rna_mode = (info_container_.mode() == config::pipeline_type::rna);

        bool use_restricted = gp_.count<SmartEdgeSet<std::unordered_set<EdgeId>, Graph>>();
        if (use_restricted) {
            DEBUG("RestrictedEdgeSet is present");
        }
        typedef std::function<bool(EdgeId edge, const std::vector<EdgeId>& path)> BulgeCallbackF;
        BulgeCallbackF bulge_callback_f = nullptr;
        if (use_restricted) {
            bulge_callback_f = [this](EdgeId e, const std::vector<EdgeId>&) {
                                   return gp_.get<SmartEdgeSet<std::unordered_set<EdgeId>, Graph>>().count(e) > 0;
                               };
        }

        INFO("Graph simplification started");
        printer_(info_printer_pos::before_simplification);

        size_t iteration = 0;
        auto message_callback = [&] () {
            INFO("PROCEDURE == Simplification cycle, iteration " << ++iteration);
        };

        CompositeAlgorithm<Graph> algo(g_, message_callback);
        auto algo_tc_br = std::make_shared<CompositeAlgorithm<Graph>>(g_);
        algo_tc_br->AddAlgo(TipClipperInstance(g_, simplif_cfg_.tc, info_container_, removal_handler_),
                            "Tip clipper");
        algo_tc_br->AddAlgo(DeadEndInstance(g_, simplif_cfg_.dead_end, info_container_, removal_handler_),
                            "Dead end clipper");
        algo_tc_br->AddAlgo(BRInstance(g_, simplif_cfg_.br, info_container_, bulge_callback_f, removal_handler_),
                            "Bulge remover");

//        algo.AddAlgo(
//                RelativelyLowCoverageDisconnectorInstance(gp_.g, gp_.flanking_cov,
//                                                          simplif_cfg_.red, info_container_),
//                "Disconnecting edges with relatively low coverage");

        if (rna_mode) {
            algo.AddAlgo<LoopedAlgorithm<Graph>>("", g_, algo_tc_br);
        } else {
            algo.AddAlgo(algo_tc_br, "");
        }

        algo.AddAlgo(ECRemoverInstance(g_, simplif_cfg_.ec, info_container_,
                                       removal_handler_),
                     "Low coverage edge remover");

        //all primary option set to closely mimic previous rna behavior
        AlgorithmRunningHelper<Graph>::IterativeThresholdsRun(algo,
                                                              simplif_cfg_.cycle_iter_count,
                                                              /*all_primary*/rna_mode);

        AlgorithmRunningHelper<Graph>::LoopedRun(algo, /*min it count*/1, /*max it count*/10);
    }
};

std::shared_ptr<visualization::graph_colorer::GraphColorer<Graph>> DefaultGPColorer(
        const GraphPack &gp) {
    auto mapper = MapperInstance(gp);
    const auto &genome = gp.get<GenomeStorage>();
    auto path1 = mapper->MapSequence(genome.GetSequence()).path();
    auto path2 = mapper->MapSequence(!genome.GetSequence()).path();
    return visualization::graph_colorer::DefaultColorer(gp.get<Graph>(), path1, path2);
}

void RawSimplification::run(GraphPack &gp, const char*) {
    using namespace omnigraph;

    //no other handlers here, todo change with DetachAll
    auto &index = gp.get_mutable<EdgeIndex<Graph>>();
    if (index.IsAttached())
        index.Detach();
    index.clear();

    visualization::graph_labeler::DefaultLabeler<Graph> labeler(gp.get<Graph>(), gp.get<EdgesPositionHandler<Graph>>());
    stats::detail_info_printer printer(gp, labeler, cfg::get().output_dir);

    GraphSimplifier simplifier(gp, CreateInfoContainer(gp),
                               preliminary_ ? *cfg::get().preliminary_simp : cfg::get().simp,
                               nullptr/*removal_handler_f*/,
                               printer);
    if (cfg::get().mode == config::pipeline_type::rna)
        simplifier.InitialCleaningRNA();
    else if (cfg::get().mode == config::pipeline_type::rnaviral) {
        simplifier.InitialCleaningRNA();
        simplifier.InitialCleaning();
    } else {
        simplifier.InitialCleaning();
    }
}

void Simplification::run(GraphPack &gp, const char*) {
    using namespace omnigraph;

    //no other handlers here, todo change with DetachAll
    auto &index = gp.get_mutable<EdgeIndex<Graph>>();
    if (index.IsAttached())
        index.Detach();
    index.clear();

    visualization::graph_labeler::DefaultLabeler<Graph> labeler(gp.get<Graph>(), gp.get<EdgesPositionHandler<Graph>>());
    stats::detail_info_printer printer(gp, labeler, cfg::get().output_dir);

    //  QualityLoggingRemovalHandler<Graph> qual_removal_handler(gp.g, edge_qual);
//    auto colorer = debruijn_graph::DefaultGPColorer(gp);
//    QualityEdgeLocalityPrintingRH<Graph> qual_removal_handler(gp.g, gp.edge_qual, labeler, colorer,
//                                   cfg::get().output_dir + "pictures/colored_edges_deleted/");
//
//    //positive quality edges removed (folder colored_edges_deleted)
//    std::function<void(EdgeId)> removal_handler_f = boost::bind(
//            //            &QualityLoggingRemovalHandler<Graph>::HandleDelete,
//            &QualityEdgeLocalityPrintingRH<Graph>::HandleDelete,
//            boost::ref(qual_removal_handler), _1);

    SimplifInfoContainer info_container = CreateInfoContainer(gp);

    GraphSimplifier simplifier(gp, info_container,
                               preliminary_ ? *cfg::get().preliminary_simp : cfg::get().simp,
                               nullptr/*removal_handler_f*/,
                               printer);
    simplifier.SimplifyGraph();
    CompressAllVertices(gp.get_mutable<Graph>());
}

void SimplificationCleanup::run(GraphPack &gp, const char*) {
    const auto &graph = gp.get<Graph>();
    visualization::graph_labeler::DefaultLabeler<Graph> labeler(graph, gp.get<EdgesPositionHandler<Graph>>());
    stats::detail_info_printer printer(gp, labeler, cfg::get().output_dir);

    GraphSimplifier simplifier(gp, CreateInfoContainer(gp),
                               cfg::get().simp,
                               nullptr/*removal_handler_f*/,
                               printer);

    simplifier.PostSimplification();

    DEBUG("Graph simplification finished");

    INFO("Counting average coverage");
    AvgCovereageCounter<Graph> cov_counter(graph);

    VERIFY(cfg::get().ds.average_coverage == 0.);
    cfg::get_writable().ds.average_coverage = cov_counter.Count();

    INFO("Average coverage = " << cfg::get().ds.average_coverage);
    if (!cfg::get().uneven_depth) {
        if (math::ls(cfg::get().ds.average_coverage, gp.get<GenomicInfo>().ec_bound()))
            WARN("The determined erroneous connection coverage threshold may be determined improperly\n");
    }
}

} //debruijn_graph
