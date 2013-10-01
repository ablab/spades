//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "standard.hpp"

#include "logger/logger.hpp"
#include "debruijn_stats.hpp"
#include "omni_labelers.hpp"

#include "de/distance_estimation.hpp"
#include "de/smoothing_distance_estimation.hpp"

#include "omni/omni_utils.hpp"

#include "path_extend/path_extend_launch.hpp"
#include "contig_output.hpp"

#include "pac_index.hpp"
#include "long_read_storage.hpp"
#include "loop_filter.hpp"
#include "coverage_based_rr.hpp"
#include "pacbio_aligner.hpp"

#include "repeat_resolving.hpp"

namespace debruijn_graph {

#if 0
void resolve_repeats_by_coverage(conj_graph_pack& conj_gp, size_t insert_size, std::vector< PathInfo<Graph> >& filteredPaths,
                                 const EdgeQuality<Graph, Index>& quality_labeler) {

    typedef DeBruijnEdgeIndex<KmerStoringDeBruijnEdgeIndex<conj_graph_pack::graph_t, runtime_k::RtSeq>> KmerIndex;

    KmerIndex kmer_index((unsigned) conj_gp.g.k() + 1, conj_gp.g, cfg::get().output_dir);
    if (cfg::get().developer_mode) {

        std::string path;
        if (cfg::get().entry_point < ws_repeats_resolving)
            path = cfg::get().output_dir + "/saves/debruijn_kmer_index_after_construction";
        else
            path = cfg::get().load_from + "/debruijn_kmer_index_after_construction";
        bool val = LoadEdgeIndex(path, kmer_index);
        VERIFY_MSG(val, "can not open file "+path+".kmidx");
        //		INFO("Updating index from graph started");
        //        EdgeInfoUpdater<KmerIndex, Graph> updater(conj_gp.g, kmer_index);
        //        updater.UpdateAll();
    }

    FlankingCoverage<conj_graph_pack::graph_t, conj_graph_pack::index_t> index(conj_gp.g, kmer_index, 50);
    EdgeLabelHandler<conj_graph_pack::graph_t> labels_after(conj_gp.g, conj_gp.g);
    auto cov_rr = CoverageBasedResolution<conj_graph_pack> (&conj_gp, cfg::get().cbrr.coverage_threshold_one_list, cfg::get().cbrr.coverage_threshold_match,
                                                            cfg::get().cbrr.coverage_threshold_global, cfg::get().cbrr.tandem_ratio_lower_threshold, cfg::get().cbrr.tandem_ratio_upper_threshold, cfg::get().cbrr.repeat_length_upper_threshold);
    cov_rr.resolve_repeats_by_coverage(index, insert_size, labels_after, quality_labeler, conj_gp.clustered_indices[0], filteredPaths);

    INFO("Repeats are resolved by coverage");
}
#endif

size_t get_first_pe_lib_index() {
    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i)
        if (cfg::get().ds.reads[i].type() == io::LibraryType::PairedEnd && cfg::get().ds.reads[i].data().mean_insert_size != 0.0)
            return i;

    return -1UL;
}

void pe_resolving(conj_graph_pack& conj_gp, const EdgeQuality<Graph, Index>& quality_labeler) {
    std::vector<PairedIndexT*> pe_indexs;
    std::vector<PairedIndexT*> pe_scaf_indexs;
    std::vector<size_t> indexes;

    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
        const auto& lib = cfg::get().ds.reads[i];
        if (lib.data().mean_insert_size != 0.0 &&
            (lib.type() == io::LibraryType::PairedEnd ||
             lib.type() == io::LibraryType::MatePairs)) {
            pe_indexs.push_back(&conj_gp.clustered_indices[i]);
            pe_scaf_indexs.push_back(&conj_gp.scaffolding_indices[i]);
            indexes.push_back(i);
        }
    }

    //LongReadStorage<Graph> long_read(conj_gp.g);
    //long_read.LoadFromFile("/storage/labnas/students/igorbunova/path-extend2/algorithmic-biology/assembler/pacbio.mpr");

    PathStorage<Graph> long_read(conj_gp.g);
    GapStorage<Graph> gaps(conj_gp.g);

    std::vector< PathInfo<Graph> > filteredPaths;
#if 0
    if (cfg::get().coverage_based_rr_on) {
        size_t pe_lib_index = get_first_pe_lib_index();
        const io::SequencingLibrary<debruijn_config::DataSetData> &lib = cfg::get().ds.reads[pe_lib_index];
        resolve_repeats_by_coverage(conj_gp, (size_t) lib.data().mean_insert_size, filteredPaths, quality_labeler);
    }
#endif

    //LongReadStorage<Graph> long_read(conj_gp.g);
    if (cfg::get().pacbio_test_on) {
        INFO("creating multiindex with k = " << cfg::get().pb.pacbio_k);
        PacBioAligner pac_aligner(conj_gp, cfg::get().pb.pacbio_k);
        INFO("index created");
        filteredPaths = long_read.GetAllPaths();
        pac_aligner.pacbio_test(long_read, gaps);
    }

    if (cfg::get().use_scaffolder && cfg::get().pe_params.param_set.scaffolder_options.on) {
        //path_extend::resolve_repeats_pe(conj_gp, pe_indexs, pe_scaf_indexs, indexes, long_read.GetAllPaths(), cfg::get().output_dir, "scaffolds.fasta", true, boost::optional<std::string>("final_contigs.fasta"));
        path_extend::resolve_repeats_pe(conj_gp, pe_indexs, pe_scaf_indexs, indexes, filteredPaths, cfg::get().output_dir, "scaffolds.fasta", true, boost::optional<std::string>("final_contigs.fasta"));
    } else {
        pe_scaf_indexs.clear();
        //path_extend::resolve_repeats_pe(conj_gp, pe_indexs, pe_scaf_indexs, indexes, long_read.GetAllPaths(), cfg::get().output_dir, "final_contigs.fasta", false, boost::none);
        path_extend::resolve_repeats_pe(conj_gp, pe_indexs, pe_scaf_indexs, indexes, filteredPaths, cfg::get().output_dir, "final_contigs.fasta", false, boost::none);
    }
}

void RepeatResolution::run(conj_graph_pack &gp) {
    OutputContigs(gp.g, cfg::get().additional_contigs, cfg::get().use_unipaths,
                  cfg::get().simp.tec.plausibility_length);
    OutputContigs(gp.g, cfg::get().output_dir + "before_rr.fasta");

    if (cfg::get().developer_mode && cfg::get().pos.late_threading) {
        FillPos(gp, gp.genome, "10");
        FillPos(gp, !gp.genome, "11");
        if (!cfg::get().pos.contigs_for_threading.empty() &&
            FileExists(cfg::get().pos.contigs_for_threading))
          FillPosWithRC(gp, cfg::get().pos.contigs_for_threading, "thr_");

        if (!cfg::get().pos.contigs_to_analyze.empty() &&
            FileExists(cfg::get().pos.contigs_to_analyze))
          FillPosWithRC(gp, cfg::get().pos.contigs_to_analyze, "anlz_");
    }

    //todo refactor labeler creation
    total_labeler_graph_struct graph_struct(gp.g, &gp.int_ids,
                                            &gp.edge_pos);
    total_labeler tot_lab(&graph_struct);
    EdgeQuality<Graph, Index> quality_labeler(gp.g, gp.index,
                                              gp.kmer_mapper, gp.genome);
    CompositeLabeler<Graph> labeler(tot_lab, quality_labeler);
    detail_info_printer printer(gp, labeler, cfg::get().output_dir);
    printer(ipp_before_repeat_resolution);

    bool no_valid_libs = true;
    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i)
        if (cfg::get().ds.reads[i].data().mean_insert_size != 0.0) {
            no_valid_libs = false;
            break;
        }

    if (cfg::get().paired_mode && no_valid_libs)
        WARN("Insert size was not estimated for any of the paired libraries, repeat resolution module will not run.");

    if (no_valid_libs ||
        cfg::get().rm == debruijn_graph::resolving_mode::rm_none) {
        OutputContigs(gp.g, cfg::get().output_dir + "final_contigs.fasta");

        if (cfg::get().pacbio_test_on) {
            PathStorage<Graph> long_read(gp.g);
            GapStorage<Graph> gaps(gp.g);
            std::vector< PathInfo<Graph> > filteredPaths;
            //LongReadStorage<Graph> long_read(gp.g);
            INFO("creating  multiindex with k = " << cfg::get().pb.pacbio_k);
            PacBioAligner pac_aligner(gp, cfg::get().pb.pacbio_k);
            INFO("index created");
            filteredPaths = long_read.GetAllPaths();
            pac_aligner.pacbio_test(long_read, gaps);
        }
        return;
    }

    // Repeat resolving begins
    size_t pe_lib_index = get_first_pe_lib_index();
    if (cfg::get().ds.reads.lib_count() > 1 || pe_lib_index == -1UL ||
        cfg::get().rm == debruijn_graph::resolving_mode::rm_path_extend) {
        INFO("Path-Extend repeat resolving");
        pe_resolving(gp, quality_labeler);
    } else {
        INFO("Unsupported repeat resolver");
        OutputContigs(gp.g, cfg::get().output_dir + "final_contigs.fasta");
    }
}

void ContigOutput::run(conj_graph_pack &gp) {
    OutputContigs(gp.g, cfg::get().additional_contigs, cfg::get().use_unipaths,
                  cfg::get().simp.tec.plausibility_length);
    OutputContigs(gp.g, cfg::get().output_dir + "before_rr.fasta");
    OutputContigs(gp.g, cfg::get().output_dir + "final_contigs.fasta");
}


} // debruijn_graph
