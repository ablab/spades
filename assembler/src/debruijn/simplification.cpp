//****************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
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
    
    if (!cfg::get().developer_mode && gp.index.IsAttached())
        gp.index.Detach();
    
    boost::function<void(EdgeId)> removal_handler_f (0);

    if (cfg::get().developer_mode) {
        CollectPositions(gp);
        gp.ClearQuality();
        gp.FillQuality();
        
        QualityLoggingRemovalHandler<Graph> qual_removal_handler(gp.g, gp.edge_qual);

//        auto colorer = debruijn_graph::DefaultGPColorer(gp);
//        QualityEdgeLocalityPrintingRH<Graph> qual_removal_handler(gp.g, gp.edge_qual, labeler, colorer,
//                                       cfg::get().output_dir + "pictures/colored_edges_deleted/");
//
//        //positive quality edges removed (folder colored_edges_deleted)
        removal_handler_f = //0
            boost::bind(
            &QualityLoggingRemovalHandler<Graph>::HandleDelete,
//                &QualityEdgeLocalityPrintingRH<Graph>::HandleDelete,
                boost::ref(qual_removal_handler), _1);

    }

    omnigraph::DefaultLabeler<Graph> labeler(gp.g, gp.edge_pos);
    stats::detail_info_printer printer(gp, labeler, cfg::get().output_dir);

    debruijn::simplification::SimplifyGraph(gp, removal_handler_f,
                  printer, /*iteration count*/10);

    AvgCovereageCounter<Graph> cov_counter(gp.g);
    cfg::get_writable().ds.set_avg_coverage(cov_counter.Count());
}

void SimplificationCleanup::run(conj_graph_pack &gp, const char*) {
    omnigraph::DefaultLabeler<Graph> labeler/*tot_lab*/(gp.g, gp.edge_pos);
    stats::detail_info_printer printer(gp, labeler, cfg::get().output_dir);

    printer(ipp_removing_isolated_edges);

    debruijn::simplification::RemoveIsolatedEdges(gp.g, cfg::get().simp.ier, cfg::get().ds.RL());
//todo return this functionality
//        INFO("Removed " << removed << " edges");

    size_t low_threshold = gp.ginfo.trusted_bound();
    if (low_threshold) {
      EdgeRemover<Graph> remover(gp.g);
      INFO("Removing all the edges having coverage " << low_threshold << " and less");
      size_t cnt = 0;
      for (auto it = gp.g.SmartEdgeBegin(); !it.IsEnd(); ++it)
        if (math::le(gp.g.coverage(*it), (double)low_threshold)) {
          remover.DeleteEdge(*it);
          cnt += 1;
        }
      INFO("Deleted " << cnt << " edges");
    }

    printer(ipp_final_simplified);

    // FIXME: Get rid of this
    if (cfg::get().correct_mismatches || cfg::get().rr_enable) {
        INFO("Final index refill");
        gp.index.Refill();
        INFO("Final index refill finished");
        if (!gp.index.IsAttached())
            gp.index.Attach();
    }

    DEBUG("Graph simplification finished");

    INFO("Counting average coverage");
    AvgCovereageCounter<Graph> cov_counter(gp.g);
    cfg::get_writable().ds.set_avg_coverage(cov_counter.Count());
    INFO("Average coverage = " << cfg::get().ds.avg_coverage());
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

void correct_mismatches(conj_graph_pack &gp) {
    INFO("Correcting mismatches");
    auto_ptr<io::IReader<io::SingleReadSeq>> stream = single_binary_multireader(true, true);
    size_t corrected = MismatchShallNotPass<conj_graph_pack, io::SingleReadSeq>(gp, 2).StopAllMismatches(*stream, 1);
    INFO("Corrected " << corrected << " nucleotides");
}

void parallel_correct_mismatches(conj_graph_pack &gp) {
    INFO("Correcting mismatches");
    auto streams = single_binary_readers(true,  true);
    size_t corrected = MismatchShallNotPass<conj_graph_pack, io::SingleReadSeq>(gp, 2).ParallelStopAllMismatches(*streams, 1);
    INFO("Corrected " << corrected << " nucleotides");
}

void exec_simplification(conj_graph_pack& gp) {
    simplify_graph(gp);

    if (cfg::get().correct_mismatches)
    {
        parallel_correct_mismatches(gp);
        }
    save_simplification(gp);
    if (cfg::get().graph_read_corr.enable) {
        //			corrected_and_save_reads(gp);
    }

    } else {
        INFO("Loading Simplification");

        path::files_t used_files;
        load_simplification(gp, &used_files);
        link_files_by_prefix(used_files, cfg::get().output_saves);
        //		if (cfg::get().correct_mismatches) {
        //			parallel_correct_mismatches(gp);
        //		}
    }
}
#endif

} //debruijn_graph
