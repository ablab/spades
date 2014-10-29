//****************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "standard.hpp"
#include "graph_simplification.hpp"
#include "omni/visualization/graph_labeler.hpp"
#include "io/single_read.hpp"
#include "positions.hpp"

#include "simplification.hpp"

namespace debruijn_graph {

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
//    boost::function<void(EdgeId)> removal_handler_f = boost::bind(
//            //            &QualityLoggingRemovalHandler<Graph>::HandleDelete,
//            &QualityEdgeLocalityPrintingRH<Graph>::HandleDelete,
//            boost::ref(qual_removal_handler), _1);
    debruijn::simplification::SimplifyGraph(gp, 0/*removal_handler_f*/,
                  printer, /*iteration count*/10);


    AvgCovereageCounter<Graph> cov_counter(gp.g);
    cfg::get_writable().ds.set_avg_coverage(cov_counter.Count());
}

void SimplificationCleanup::run(conj_graph_pack &gp, const char*) {
    omnigraph::DefaultLabeler<Graph> labeler/*tot_lab*/(gp.g, gp.edge_pos);
    stats::detail_info_printer printer(gp, labeler, cfg::get().output_dir);

    printer(ipp_removing_isolated_edges);

    debruijn::simplification::RemoveIsolatedEdges(gp.g, cfg::get().simp.ier, cfg::get().ds.RL(), boost::function<void(EdgeId)>(0), cfg::get().max_threads);

    double low_threshold = gp.ginfo.trusted_bound();
    if (math::ge(low_threshold, 0.0)) {
        INFO("Removing all the edges having coverage " << low_threshold << " and less");
        omnigraph::EdgeRemovingAlgorithm<Graph> removing_algo(gp.g,
                                                              std::make_shared<func::AlwaysTrue<EdgeId>>(), 0);

        removing_algo.Process(CoverageComparator<Graph>(gp.g),
                              std::make_shared<CoverageUpperBound<Graph>>(gp.g, low_threshold));
    }

    printer(ipp_final_simplified);

    DEBUG("Graph simplification finished");

    INFO("Counting average coverage");
    AvgCovereageCounter<Graph> cov_counter(gp.g);
    cfg::get_writable().ds.set_avg_coverage(cov_counter.Count());
    INFO("Average coverage = " << cfg::get().ds.avg_coverage());
    if (!cfg::get().ds.single_cell) {
        if (cfg::get().ds.avg_coverage() < gp.ginfo.ec_bound())
            WARN("The determined errorneous connection coverage threshold may be dtermined improperly\n");
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
