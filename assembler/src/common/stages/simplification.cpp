//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "assembly_graph/core/basic_graph_stats.hpp"
#include "assembly_graph/graph_support/graph_processing_algorithm.hpp"
#include "stages/simplification_pipeline/simplification_settings.hpp"
#include "stages/simplification_pipeline/graph_simplification.hpp"
#include "modules/simplification/parallel_simplification_algorithms.hpp"

#include "simplification.hpp"

namespace debruijn_graph {

using namespace debruijn::simplification;
using namespace config;

template<class graph_pack>
shared_ptr<visualization::graph_colorer::GraphColorer<typename graph_pack::graph_t>> DefaultGPColorer(
    const graph_pack& gp) {
    io::SingleRead genome("ref", gp.genome.str());
    auto mapper = MapperInstance(gp);
    auto path1 = mapper->MapRead(genome).path();
    auto path2 = mapper->MapRead(!genome).path();
    return visualization::graph_colorer::DefaultColorer(gp.g, path1, path2);
}

class GraphSimplifier {
    typedef std::function<void(EdgeId)> HandlerF;

    typedef std::vector<std::pair<AlgoPtr<Graph>, std::string>> AlgoStorageT;

    conj_graph_pack& gp_;
    Graph& g_;
    SimplifInfoContainer info_container_;
    const debruijn_config::simplification simplif_cfg_;

    CountingCallback<Graph> cnt_callback_;
    HandlerF removal_handler_;
    stats::detail_info_printer& printer_;

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

    void RemoveShortPolyATEdges(size_t max_length,
                                HandlerF removal_handler = 0, size_t chunk_cnt = 1) {
        INFO("Removing short polyAT");
        EdgeRemover<Graph> er(g_, removal_handler);
        ATCondition<Graph> condition (g_, 0.8, max_length, false);
        for (auto iter = g_.SmartEdgeBegin(); !iter.IsEnd(); ++iter){
            if (g_.length(*iter) == 1 && condition.Check(*iter)) {
                er.DeleteEdgeNoCompress(*iter);
            }
        }
        ParallelCompress(g_, chunk_cnt);
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

        if (info_container_.mode() == config::pipeline_type::rna){
            RemoveShortPolyATEdges(1, removal_handler_, info_container_.chunk_cnt());
            PushValid(ShortPolyATEdgesRemoverInstance(g_, 1, removal_handler_, info_container_.chunk_cnt()), "Short PolyA/T Edges",algos) ;
            PushValid(ATTipClipperInstance(g_, removal_handler_, info_container_.chunk_cnt()), "AT Tips", algos);
        }

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

            PushValid(
                    LowFlankDisconnectorInstance(g_, gp_.flanking_cov,
                                                 simplif_cfg_.init_clean.disconnect_flank_cov, info_container_,
                                                 removal_handler_),
                    "Disconnecting edges with low flanking coverage",
                    algos);
        }

        RunAlgos(algos);

        if (info_container_.mode() == config::pipeline_type::rna){
            RemoveHiddenLoopEC(g_, gp_.flanking_cov, info_container_.detected_coverage_bound(), simplif_cfg_.her, removal_handler_);
            cnt_callback_.Report();
        }
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
        bool changed = false;
        if (simplif_cfg_.topology_simplif_enabled && info_container_.main_iteration()) {
            changed |= AllTopology();
            changed |= MaxFlowRemoveErroneousEdges(gp_.g, simplif_cfg_.mfec,
                                                   removal_handler_);
            cnt_callback_.Report();
        }
        return changed;
    }

    void PostSimplification() {
        using namespace omnigraph;
        using namespace pred;
        INFO("PROCEDURE == Post simplification");
        size_t iteration = 0;

        AlgoStorageT algos;

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


        typename ComponentRemover<Graph>::HandlerF set_removal_handler_f;
        if (removal_handler_) {
            set_removal_handler_f = [=](const set<EdgeId>& edges) {
                std::for_each(edges.begin(), edges.end(), removal_handler_);
            };
        }

        PushValid(
                RelativelyLowCoverageComponentRemoverInstance(gp_.g, gp_.flanking_cov,
                                                              simplif_cfg_.rcc, info_container_, set_removal_handler_f),
                "Relative coverage component remover",
                algos);


        PushValid(
                RelativelyLowCoverageDisconnectorInstance(gp_.g, gp_.flanking_cov,
                                                          simplif_cfg_.relative_ed, info_container_),
                "Disconnecting edges with relatively low coverage",
                algos);

        PushValid(
                ComplexTipClipperInstance(gp_.g, simplif_cfg_.complex_tc, info_container_, set_removal_handler_f),
                "Complex tip clipper",
                algos);

        PushValid(
                ComplexBRInstance(gp_.g, simplif_cfg_.cbr, info_container_, iteration),
                "Complex bulge remover",
                algos);

        PushValid(
                TipClipperInstance(g_, simplif_cfg_.tc,
                                   info_container_, removal_handler_),
                "Tip clipper",
                algos);

        PushValid(
                TipClipperInstance(g_, simplif_cfg_.final_tc,
                                   info_container_, removal_handler_),
                "Final tip clipper",
                algos);

        PushValid(
                BRInstance(g_, simplif_cfg_.br,
                                   info_container_, removal_handler_),
                "Bulge remover",
                algos);

        PushValid(
                BRInstance(g_, simplif_cfg_.final_br,
                                   info_container_, removal_handler_),
                "Final bulge remover",
                algos);

        if (simplif_cfg_.topology_simplif_enabled) {
            PushValid(
                    TopologyTipClipperInstance(g_, simplif_cfg_.ttc,
                                                      info_container_, removal_handler_),
                    "Topology tip clipper",
                    algos);
        }

        //FIXME need better configuration

        if (info_container_.mode() == config::pipeline_type::meta) {
            PushValid(
                    BRInstance(g_, simplif_cfg_.second_final_br,
                                       info_container_, removal_handler_),
                    "Yet another final bulge remover",
                    algos);

            EdgePredicate<Graph> meta_thorn_condition
                    = And(LengthUpperBound<Graph>(g_, LengthThresholdFinder::MaxErroneousConnectionLength(
                                                                           g_.k(), simplif_cfg_.isec.max_ec_length_coefficient)),

                      And([&] (EdgeId e) {
                              //todo configure!
                              return simplification::relative_coverage::
                                         RelativeCoverageHelper<Graph>(g_, gp_.flanking_cov, 2).AnyHighlyCoveredOnFourSides(e);
                          },

                      And(UniqueIncomingPathLengthLowerBound(g_, simplif_cfg_.isec.uniqueness_length),

                          //todo configure!
                          TopologicalThornCondition<Graph>(g_, simplif_cfg_.isec.span_distance, /*max edge cnt*/5))));

            PushValid(std::make_shared<ParallelEdgeRemovingAlgorithm<Graph>>(g_, meta_thorn_condition, info_container_.chunk_cnt(), 
                      removal_handler_),
                      "Thorn remover (meta)",
                      algos);
        }

        if (info_container_.mode() == config::pipeline_type::rna) {
            PushValid(ATTipClipperInstance(g_, removal_handler_, info_container_.chunk_cnt()), "AT Tips", algos);
        }

        bool enable_flag = true;
        while (enable_flag) {
            enable_flag = false;

            INFO("Iteration " << iteration);

            enable_flag |= FinalRemoveErroneousEdges();
            cnt_callback_.Report();

            enable_flag |= RunAlgos(algos);

            iteration++;

            //    printer(ipp_before_final_err_con_removal);
            //        printer(ipp_final_tip_clipping, str(format("_%d") % iteration));
            //        printer(ipp_final_err_con_removal, str(format("_%d") % iteration));
            //        printer(ipp_final_bulge_removal, str(format("_%d") % iteration));
        }

        if (simplif_cfg_.topology_simplif_enabled) {
            RemoveHiddenEC(gp_.g, gp_.flanking_cov, simplif_cfg_.her, info_container_, removal_handler_);

            cnt_callback_.Report();
        }

        if (info_container_.mode() == config::pipeline_type::meta && simplif_cfg_.her.enabled) {
            VERIFY(math::ls(simplif_cfg_.her.unreliability_threshold, 0.));
            MetaHiddenECRemover<Graph> algo(g_, info_container_.chunk_cnt(), gp_.flanking_cov,
                                                      simplif_cfg_.her.uniqueness_length,
                                                      simplif_cfg_.her.relative_threshold,
                                                      removal_handler_);
            INFO("Running Hidden EC remover (meta)");
            LoopedRun(algo);
            cnt_callback_.Report();
        }

        INFO("Disrupting self-conjugate edges");
        SelfConjugateDisruptor<Graph>(gp_.g, removal_handler_).Run();
        cnt_callback_.Report();
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

    void PushValid(const AlgoPtr<Graph>& algo_ptr, std::string comment, AlgoStorageT& algos) const {
        if (algo_ptr) {
            algos.push_back(std::make_pair(algo_ptr, comment));
        }
    }

    bool RunAlgos(AlgoStorageT& algos, bool force_primary_launch = false) {
        bool changed = false;
        for (auto algo_comment : algos) {
             INFO("Running " << algo_comment.second);
             size_t triggered = algo_comment.first->Run(force_primary_launch);
             INFO("Triggered " << triggered << " times");
             changed |= (triggered > 0);
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
        printer_(info_printer_pos::before_simplification);
        INFO("Graph simplification started");

        InitialCleaning();

        AlgoStorageT algos;

        PushValid(
                TipClipperInstance(g_, simplif_cfg_.tc, info_container_, removal_handler_, simplif_cfg_.cycle_iter_count),
                "Tip clipper",
                algos);
        PushValid(
                BRInstance(g_, simplif_cfg_.br, info_container_, removal_handler_, simplif_cfg_.cycle_iter_count),
                "Bulge remover",
                algos);
        PushValid(
                ECRemoverInstance(g_, simplif_cfg_.ec, info_container_, removal_handler_, simplif_cfg_.cycle_iter_count),
                "Low coverage edge remover",
                algos);

        size_t iteration = 0;
        bool graph_changed = true;
        //cannot stop simply if nothing changed, since threshold change on every iteration
        while (iteration < simplif_cfg_.cycle_iter_count || graph_changed) {
            INFO("PROCEDURE == Simplification cycle, iteration " << iteration + 1);
            graph_changed = RunAlgos(algos);
            ++iteration;
        }

        printer_(info_printer_pos::before_post_simplification);

        if (simplif_cfg_.post_simplif_enabled) {
            PostSimplification();
        } else {
            INFO("PostSimplification disabled");
        }
    }

    //FIXME reduce code duplication
    void SimplifyRNAGraph() {
        printer_(info_printer_pos::before_simplification);
        INFO("Graph simplification started");

        InitialCleaning();

        if (gp_.genome.GetSequence().size() > 0) {
            DEBUG("Reference genome length = " + std::to_string(gp_.genome.GetSequence().size()));
        }

        AlgoStorageT ec_algo;

        PushValid(ECRemoverInstance(g_, simplif_cfg_.ec, info_container_, removal_handler_,
                                            simplif_cfg_.cycle_iter_count), "Low coverage edge remover", ec_algo);

        size_t iteration = 0;
        bool graph_changed_ec = true;
        //TODO: config. Or just graph_changed?
        size_t tc_max_iteration = 2;
        //cannot stop simply if nothing changed, since threshold change on every iteration
        while (iteration < simplif_cfg_.cycle_iter_count || graph_changed_ec) {
            AlgoStorageT algos;
            PushValid(
                    TipClipperInstance(g_, simplif_cfg_.tc, info_container_, removal_handler_, tc_max_iteration),
                    "Tip clipper",
                    algos);
            PushValid(
                    DeadEndInstance(g_, simplif_cfg_.dead_end, info_container_, removal_handler_, tc_max_iteration),
                    "Dead end clipper",
                    algos);
            PushValid(
                    BRInstance(g_, simplif_cfg_.br, info_container_, removal_handler_, tc_max_iteration),
                    "Bulge remover",
                    algos);
            bool graph_changed = true;
            size_t tc_iteration = 0;

            while (tc_iteration < tc_max_iteration || graph_changed) {
                INFO("PROCEDURE == Tip clipper and bulge removal cycle, iteration " << iteration + 1 << "." << tc_iteration);
                graph_changed = RunAlgos(algos);
                ++tc_iteration;
            }
            INFO("PROCEDURE == Erroneous connection, iteration " << iteration + 1);
            graph_changed_ec = RunAlgos(ec_algo);
            ++iteration;
        }

        printer_(info_printer_pos::before_post_simplification);

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

    visualization::graph_labeler::DefaultLabeler<Graph> labeler(gp.g, gp.edge_pos);
    
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


    SimplifInfoContainer info_container(cfg::get().mode);
    info_container.set_read_length(cfg::get().ds.RL())
        .set_main_iteration(cfg::get().main_iteration)
        .set_chunk_cnt(5 * cfg::get().max_threads);

    //0 if model didn't converge
    //todo take max with trusted_bound
    //FIXME add warning when used for uneven coverage applications
    info_container.set_detected_mean_coverage(gp.ginfo.estimated_mean())
            .set_detected_coverage_bound(gp.ginfo.ec_bound());

    GraphSimplifier simplifier(gp, info_container,
                               preliminary_ ? *cfg::get().preliminary_simp : cfg::get().simp,
                               nullptr/*removal_handler_f*/,
                               printer);
    if (cfg::get().mode == pipeline_type::rna)
        simplifier.SimplifyRNAGraph();
    else
        simplifier.SimplifyGraph();

}


void SimplificationCleanup::run(conj_graph_pack &gp, const char*) {
    SimplifInfoContainer info_container(cfg::get().mode);
    info_container
        .set_read_length(cfg::get().ds.RL())
        .set_main_iteration(cfg::get().main_iteration)
        .set_chunk_cnt(5 * cfg::get().max_threads);


    auto isolated_edge_remover =
        IsolatedEdgeRemoverInstance(gp.g, cfg::get().simp.ier, info_container, (EdgeRemovalHandlerF<Graph>)nullptr);
    if (isolated_edge_remover != nullptr)
        isolated_edge_remover->Run();

    double low_threshold = gp.ginfo.trusted_bound();
    if (math::gr(low_threshold, 0.0)) {
        INFO("Removing all the edges having coverage " << low_threshold << " and less");
        ParallelEdgeRemovingAlgorithm<Graph, CoverageComparator<Graph>>
                cov_cleaner(gp.g,
                            CoverageUpperBound<Graph>(gp.g, low_threshold),
                            info_container.chunk_cnt(),
                            (EdgeRemovalHandlerF<Graph>)nullptr,
                            /*canonical_only*/true,
                            CoverageComparator<Graph>(gp.g));
        cov_cleaner.Run();
    }

    visualization::graph_labeler::DefaultLabeler<Graph> labeler(gp.g, gp.edge_pos);
    stats::detail_info_printer printer(gp, labeler, cfg::get().output_dir);
    printer(info_printer_pos::final_simplified);

    DEBUG("Graph simplification finished");

    INFO("Counting average coverage");
    AvgCovereageCounter<Graph> cov_counter(gp.g);

    cfg::get_writable().ds.set_avg_coverage(cov_counter.Count());

    INFO("Average coverage = " << cfg::get().ds.avg_coverage());
    if (!cfg::get().uneven_depth) {
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
