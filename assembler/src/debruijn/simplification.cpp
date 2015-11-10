//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "standard.hpp"
#include "simplification/simplification_settings.hpp"
#include "simplification/graph_simplification.hpp"
#include "omni/visualization/graph_labeler.hpp"
#include "io/single_read.hpp"
#include "positions.hpp"

#include "simplification.hpp"

namespace debruijn_graph {

using namespace debruijn::simplification;

class GraphSimplifier {
    typedef std::function<void(EdgeId)> HandlerF;
    typedef omnigraph::PersistentEdgeRemovingAlgorithm<Graph,
            omnigraph::ParallelInterestingElementFinder<Graph, EdgeId>,
            LengthComparator<Graph>> TipClipperT;
    typedef omnigraph::PersistentEdgeRemovingAlgorithm<Graph,
            omnigraph::ParallelInterestingElementFinder<Graph, EdgeId>,
            CoverageComparator<Graph>> ECRemoverT;

    typedef std::vector<std::pair<AlgoPtr<Graph>, std::string>> AlgoStorageT;

    conj_graph_pack& gp_;
    Graph& g_;
    SimplifInfoContainer info_container_;
    const debruijn_config::simplification simplif_cfg_;

    CountingCallback<Graph> cnt_callback_;
    HandlerF removal_handler_;
    stats::detail_info_printer& printer_;

//    bool FastModeAvailable(const SimplifInfoContainer& info, double activation_cov_threshold) {
//        const auto& cfg = cfg::get();
//
//        //todo fix logic
//        //also handles meta case for now
//        if (cfg.ds.single_cell) {
//            return !cfg::get().main_iteration;
//        }
//
//        if (math::eq(info.detected_mean_coverage(), 0.) &&
//            !cfg.kcm.use_coverage_threshold) {
//            WARN("Mean coverage wasn't reliably estimated");
//            return false;
//        }
//
//        //todo review logic
//        if (math::ls(info.detected_mean_coverage(), activation_cov_threshold) &&
//            !(cfg.kcm.use_coverage_threshold &&
//              math::ge(cfg.kcm.coverage_threshold, activation_cov_threshold))) {
//            INFO("Estimated mean coverage " << info.detected_mean_coverage() <<
//                 " is less than fast mode activation coverage " << activation_cov_threshold);
//            return false;
//        }
//
//        return true;
//    }

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

    void InitialCleaning() {
        INFO("PROCEDURE == InitialCleaning");

        AlgoStorageT algos;

        PushValid(
                SelfConjugateEdgeRemoverInstance(g_,
                                                 simplif_cfg_.init_clean.self_conj_condition,
                                                 info_container_, removal_handler_),
                "Self conjugate edge remover",
                algos);

        if (PerformInitCleaning()) {
            PushValid(
                    IsolatedEdgeRemoverInstance(g_,
                                        simplif_cfg_.init_clean.ier,
                                        info_container_, removal_handler_),
                    "Initial isolated edge remover",
                    algos);

            PushValid(
                    TipClipperInstance(g_,
                               debruijn_config::simplification::tip_clipper(simplif_cfg_.init_clean.tip_condition),
                               info_container_,
                               removal_handler_),
                    "Initial tip clipper",
                    algos);

            PushValid(
                    ECRemoverInstance(g_,
                              debruijn_config::simplification::erroneous_connections_remover(simplif_cfg_.init_clean.ec_condition),
                              info_container_,
                              removal_handler_),
                    "Initial ec remover",
                    algos);
        }

        RunAlgos(algos);
    }

    bool AllTopology() {
        bool res = TopologyRemoveErroneousEdges(gp_.g, simplif_cfg_.tec,
                                                removal_handler_);
        cnt_callback_.Report();
        res |= TopologyReliabilityRemoveErroneousEdges(gp_.g, simplif_cfg_.trec,
                                                       removal_handler_);
        cnt_callback_.Report();
        res |= RemoveThorns(gp_.g, simplif_cfg_.isec, removal_handler_);
        cnt_callback_.Report();
        res |= MultiplicityCountingRemoveErroneousEdges(gp_.g, simplif_cfg_.tec,
                                                        removal_handler_);
        cnt_callback_.Report();
        return res;
    }

    bool FinalRemoveErroneousEdges() {

    //    gp.ClearQuality();
    //    gp.FillQuality();
    //    auto colorer = debruijn_graph::DefaultGPColorer(gp);
    //    omnigraph::DefaultLabeler<typename gp_t::graph_t> labeler(gp.g, gp.edge_pos);
    //    QualityEdgeLocalityPrintingRH<Graph> qual_removal_handler(gp.g, gp.edge_qual, labeler, colorer,
    //                                   cfg::get().output_dir + "pictures/colored_edges_deleted/");
    //
    //    //positive quality edges removed (folder colored_edges_deleted)
    //    std::function<void(EdgeId)> qual_removal_handler_f = boost::bind(
    //            //            &QualityLoggingRemovalHandler<Graph>::HandleDelete,
    //            &QualityEdgeLocalityPrintingRH<Graph>::HandleDelete,
    //            boost::ref(qual_removal_handler), _1);
    //
    //    std::function<void(set<EdgeId>)> set_removal_handler_f = boost::bind(
    //                &omnigraph::simplification::SingleEdgeAdapter<set<EdgeId>>, _1, qual_removal_handler_f);
    //

        std::function<void(set<EdgeId>)> set_removal_handler_f(0);
        if (removal_handler_) {
            set_removal_handler_f = std::bind(
                &omnigraph::simplification::SingleEdgeAdapter<set<EdgeId>>, std::placeholders::_1, removal_handler_);
        }

        bool changed = RemoveRelativelyLowCoverageComponents(gp_.g, gp_.flanking_cov,
                                              simplif_cfg_.rcc, info_container_, set_removal_handler_f);

        cnt_callback_.Report();

        changed |= DisconnectRelativelyLowCoverageEdges(gp_.g, gp_.flanking_cov, simplif_cfg_.relative_ed);

        if (simplif_cfg_.topology_simplif_enabled && cfg::get().main_iteration) {
            changed |= AllTopology();
            changed |= MaxFlowRemoveErroneousEdges(gp_.g, simplif_cfg_.mfec,
                                                   removal_handler_);
            cnt_callback_.Report();
        }
        return changed;
    }

    void PostSimplification() {
        INFO("PROCEDURE == Post simplification");
        size_t iteration = 0;

        AlgoStorageT algos;

        PushValid(
                TipClipperInstance(g_, simplif_cfg_.tc,
                                   info_container_, removal_handler_),
                "Tip clipper",
                algos);

        PushValid(
                ParallelBRInstance(g_, simplif_cfg_.br,
                                   info_container_, removal_handler_),
                "Bulge remover",
                algos);

        if (simplif_cfg_.topology_simplif_enabled) {
            PushValid(
                    TopologyTipClipperInstance(g_, simplif_cfg_.ttc,
                                                      info_container_, removal_handler_),
                    "Topology tip clipper",
                    algos);
        }

        //FIXME need better configuration
        if (cfg::get().ds.meta && info_container_.main_iteration()) {
            PushValid(
                    TipClipperInstance(g_, simplif_cfg_.final_tc,
                                       info_container_, removal_handler_),
                    "Final tip clipper",
                    algos);

            PushValid(
                    ParallelBRInstance(g_, simplif_cfg_.final_br,
                                       info_container_, removal_handler_),
                    "Final bulge remover",
                    algos);

            PushValid(
                    ParallelBRInstance(g_, simplif_cfg_.second_final_br,
                                       info_container_, removal_handler_),
                    "Yet another final bulge remover",
                    algos);
        }

        bool enable_flag = true;
        while (enable_flag) {
            enable_flag = false;

            INFO("Iteration " << iteration);

            enable_flag |= FinalRemoveErroneousEdges();
            cnt_callback_.Report();

            enable_flag |=  ClipComplexTips(gp_.g, simplif_cfg_.complex_tc, removal_handler_);
            cnt_callback_.Report();

            enable_flag |= RemoveComplexBulges(gp_.g, simplif_cfg_.cbr, iteration);
            cnt_callback_.Report();

            enable_flag |= RunAlgos(algos);

            iteration++;

            //    printer(ipp_before_final_err_con_removal);
            //        printer(ipp_final_tip_clipping, str(format("_%d") % iteration));
            //        printer(ipp_final_err_con_removal, str(format("_%d") % iteration));
            //        printer(ipp_final_bulge_removal, str(format("_%d") % iteration));
        }

        //fixme move to AllTopology?
        if (simplif_cfg_.topology_simplif_enabled) {
            RemoveHiddenEC(gp_.g, gp_.flanking_cov, simplif_cfg_.her, info_container_, removal_handler_);

            cnt_callback_.Report();
        }
    }

    //inline
    //void IdealSimplification(Graph& graph,
    //                         std::function<double(EdgeId)> quality_handler_f) {
    //    for (auto iterator = graph.SmartEdgeBegin(); !iterator.IsEnd();
    //         ++iterator) {
    //        if (math::eq(quality_handler_f(*iterator), 0.))
    //            graph.DeleteEdge(*iterator);
    //    }
    //    CompressAllVertices(graph);
    //}

//    std::shared_ptr<Predicate<EdgeId>> ParseCondition(const string& condition) const {
//        ConditionParser<Graph> parser(g_, condition, info_container_);
//        return parser();
//    }

    void PushValid(const AlgoPtr<Graph>& algo_ptr, std::string comment, AlgoStorageT& algos) const {
        if (algo_ptr) {
            algos.push_back(std::make_pair(algo_ptr, comment));
        }
    }

    bool RunAlgos(AlgoStorageT& algos, bool force_primary_launch = false) {
        bool changed = false;
        for (auto algo_comment : algos) {
             INFO("Running " << algo_comment.second);
             changed |= algo_comment.first->Run(force_primary_launch);
             cnt_callback_.Report();
         }
        return changed;
    }

public:
    GraphSimplifier(conj_graph_pack &gp, const SimplifInfoContainer& info_container,
                    const debruijn_config::simplification& simplif_cfg,
                    const std::function<void(EdgeId)>& removal_handler,
                    stats::detail_info_printer& printer)
            : gp_(gp),
              g_(gp_.g),
              info_container_(info_container),
              simplif_cfg_(simplif_cfg),
              removal_handler_(AddCountingCallback(cnt_callback_, removal_handler)),
              printer_(printer) {

    }

    void SimplifyGraph() {
        printer_(ipp_before_simplification);
        INFO("Graph simplification started");

        InitialCleaning();

        AlgoStorageT algos;

        PushValid(
                TipClipperInstance(g_, simplif_cfg_.tc, info_container_, removal_handler_, simplif_cfg_.cycle_iter_count),
                "Tip clipper",
                algos);
        PushValid(
                ParallelBRInstance(g_, simplif_cfg_.br, info_container_, removal_handler_, simplif_cfg_.cycle_iter_count),
                "Bulge remover",
                algos);
        PushValid(
                ECRemoverInstance(g_, simplif_cfg_.ec, info_container_, removal_handler_, simplif_cfg_.cycle_iter_count),
                "Low coverage edge remover",
                algos);

        for (size_t i = 0; i < simplif_cfg_.cycle_iter_count; i++) {
            INFO("PROCEDURE == Simplification cycle, iteration " << i + 1);
            RunAlgos(algos);
            //cannot stop even if nothing changed, since threshold change on every iteration
        }

        printer_(ipp_before_post_simplification);

        if (simplif_cfg_.post_simplif_enabled) {
            PostSimplification();
        } else {
            INFO("PostSimplification disabled");
        }
    }
};

void Simplification::run(conj_graph_pack &gp, const char*) {
    using namespace omnigraph;

    //no other handlers here, todo change with DetachAll
    gp.index.Detach();
    gp.index.clear();

    omnigraph::DefaultLabeler<Graph> labeler(gp.g, gp.edge_pos);
    
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


    SimplifInfoContainer info_container;
    info_container.set_read_length(cfg::get().ds.RL())
        .set_main_iteration(cfg::get().main_iteration)
        .set_chunk_cnt(5 * cfg::get().max_threads);

    if (!cfg::get().ds.meta) {
        //0 if model didn't converge
        //todo take max with trusted_bound
        info_container.set_detected_mean_coverage(gp.ginfo.estimated_mean())
                .set_detected_coverage_bound(gp.ginfo.ec_bound());
    }

    debruijn::simplification::GraphSimplifier simplifier(gp, info_container,
                                                                 preliminary_ ? cfg::get().preliminary_simp : cfg::get().simp,
                                                                 nullptr/*removal_handler_f*/,
                                                                 printer);
    simplifier.SimplifyGraph();
}

void SimplificationCleanup::run(conj_graph_pack &gp, const char*) {
    omnigraph::DefaultLabeler<Graph> labeler/*tot_lab*/(gp.g, gp.edge_pos);
    stats::detail_info_printer printer(gp, labeler, cfg::get().output_dir);

    SimplifInfoContainer info_container;
    info_container
        .set_read_length(cfg::get().ds.RL())
        .set_main_iteration(cfg::get().main_iteration)
        .set_chunk_cnt(5 * cfg::get().max_threads);

    IsolatedEdgeRemoverInstance(gp.g, cfg::get().simp.ier, info_container, (HandlerF<Graph>)nullptr)->Run();

    double low_threshold = gp.ginfo.trusted_bound();
    if (math::gr(low_threshold, 0.0)) {
        INFO("Removing all the edges having coverage " << low_threshold << " and less");

        ParallelEdgeRemovingAlgorithm<Graph, CoverageComparator<Graph>> cov_cleaner(gp.g,
                                                std::make_shared<CoverageUpperBound<Graph>>(gp.g, low_threshold),
                                                info_container.chunk_cnt(),
                                                (HandlerF<Graph>)nullptr,
                                                /*canonical_only*/true,
                                                CoverageComparator<Graph>(gp.g));
        cov_cleaner.Run();
    }

    printer(ipp_final_simplified);

    DEBUG("Graph simplification finished");

    INFO("Counting average coverage");
    AvgCovereageCounter<Graph> cov_counter(gp.g);
    cfg::get_writable().ds.set_avg_coverage(cov_counter.Count());
    INFO("Average coverage = " << cfg::get().ds.avg_coverage());
    if (!cfg::get().ds.single_cell) {
        if (cfg::get().ds.avg_coverage() < gp.ginfo.ec_bound())
            WARN("The determined erroneous connection coverage threshold may be determined improperly\n");
    }
}


#if 0
void corrected_and_save_reads(const conj_graph_pack& gp) {
    //saving corrected reads
    //todo read input files, correct, save and use on the next iteration

    auto_ptr<io::IReader<io::PairedReadSeq>> paired_stream =
            paired_binary_multireader(false, /*insert_size*/0);
    io::ModifyingWrapper<io::PairedReadSeq> refined_paired_stream(
        *paired_stream,
        GraphReadCorrectorInstance(gp.g, *MapperInstance(gp)));

    auto_ptr<io::IReader<io::SingleReadSeq>> single_stream =
            single_binary_multireader(false, /*include_paired_reads*/false);
    io::ModifyingWrapper<io::SingleReadSeq> refined_single_stream(
        *single_stream,
        GraphReadCorrectorInstance(gp.g, *MapperInstance(gp)));

    if (cfg::get().graph_read_corr.binary) {
        INFO("Correcting paired reads");

        io::BinaryWriter paired_converter(
            cfg::get().paired_read_prefix + "_cor", cfg::get().max_threads,
            cfg::get().buffer_size);
        paired_converter.ToBinary(refined_paired_stream);

        INFO("Correcting single reads");
        io::BinaryWriter single_converter(
            cfg::get().single_read_prefix + "_cor", cfg::get().max_threads,
            cfg::get().buffer_size);
        single_converter.ToBinary(refined_single_stream);
    } else {
        //save in fasta
        VERIFY(false);
    }

    INFO("Error correction done");
}
#endif

} //debruijn_graph
