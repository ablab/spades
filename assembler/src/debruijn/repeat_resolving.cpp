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

#include "long_read_storage.hpp"

#if 0
#include "loop_filter.hpp"
#include "pac_index.hpp"
#include "coverage_based_rr.hpp"
#include "pacbio_aligner.hpp"
#endif

#include "repeat_resolving.hpp"

namespace debruijn_graph {

#if 0
void resolve_repeats_by_coverage(conj_graph_pack& conj_gp, size_t insert_size, std::vector< PathInfo<Graph> >& filteredPaths,
                                 PairedIndexT& clustered_index,
                                 const EdgeQuality<Graph, Index>& quality_labeler) {

    typedef DeBruijnEdgeIndex<KmerStoringDeBruijnEdgeIndex<conj_graph_pack::graph_t, runtime_k::RtSeq>> KmerIndex;
        KmerIndex kmer_index((unsigned) conj_gp.g.k() + 1, conj_gp.g, cfg::get().output_dir);
        if (cfg::get().developer_mode) {

        std::string path;
        if (cfg::get().entry_point <= ws_simplification)
            path = cfg::get().output_dir + "/saves/debruijn_kmer_index_after_construction";
        else
            path = cfg::get().load_from + "/debruijn_kmer_index_after_construction";
        bool val = LoadEdgeIndex(path, kmer_index);
        VERIFY_MSG(val, "can not open file "+path+".kmidx");
        INFO("Updating index from graph");

        //DeBruijnEdgeIndexBuilder<runtime_k::RtSeq>().UpdateIndexFromGraph(kmer_index, conj_gp.g);
        EdgeInfoUpdater<KmerIndex, Graph> updater(conj_gp.g, kmer_index);
        updater.UpdateAll();
        SaveEdgeIndex(cfg::get().output_dir + "/saves/debruijn_kmer_index_after_construction", kmer_index);
    }

    /*int number_of_buckets = 10;
    auto bm = BucketMapper<conj_graph_pack::graph_t>(conj_gp.g, kmerIndex, cfg::get().K + 1, number_of_buckets);
    bm.InitBuckets();

    int bucket_in = 0;
    int bucket_out = 0;
    int repeat_distance = 500;
    double probability = bm.GetProbablityFromBucketToBucketForDistance (bucket_in, bucket_out, repeat_distance) ;*/
    //auto index = FlankingCoverage<Graph>(conj_gp.g, kmer_index, 50, cfg::get().K + 1);
    FlankingCoverage<Graph, KmerIndex> index(conj_gp.g, kmer_index, 50);
    EdgeLabelHandler<conj_graph_pack::graph_t> labels_after(conj_gp.g, conj_gp.g);
    CoverageBasedResolution<conj_graph_pack, const EdgeQuality<Graph, Index>, KmerIndex> cov_rr
    (conj_gp, kmer_index, quality_labeler, cfg::get().cbrr.tandem_ratio_lower_threshold,
            cfg::get().cbrr.tandem_ratio_upper_threshold, cfg::get().cbrr.repeat_length_upper_threshold);
    cov_rr.resolve_repeats_by_coverage(index, insert_size, labels_after, clustered_index, filteredPaths, cfg::get().output_dir + "resolved_by_coverage.fasta");

    INFO("Repeats are resolved by coverage");

}
#endif


/*
 * Return index of first paired-end library or -1 if there is no paired end library
 */
size_t GetFirstPELibIndex() {
    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i)
        if (cfg::get().ds.reads[i].type() == io::LibraryType::PairedEnd
                && cfg::get().ds.reads[i].data().mean_insert_size != 0.0)
            return i;
    return -1UL;
}

//TODO: get rid of this conversion
void ConvertLongReads(LongReadContainerT& single_long_reads, vector<PathStorageInfo<Graph> > &long_reads_libs) {
    for (size_t i = 0; i < single_long_reads.size(); ++i) {
        INFO("converting " << i)
        PathStorage<Graph>& storage = single_long_reads[i];
        vector<PathInfo<Graph> > paths = storage.GetAllPaths();
        PathStorageInfo<Graph> single_storage(paths,
                cfg::get().pe_params.long_reads.single_reads.filtering,
                cfg::get().pe_params.long_reads.single_reads.weight_priority,
                cfg::get().pe_params.long_reads.single_reads.unique_edge_priority);
        long_reads_libs.push_back(single_storage);
        INFO("done " << i)
    }
}

void pe_resolving(conj_graph_pack& gp, const EdgeQuality<Graph, Index>& /* quality_labeler */) {
    vector<PairedIndexT*> pe_indexes;
    vector<PairedIndexT*> pe_scaf_indices;
    vector<size_t> indexes;
    vector<PathStorageInfo<Graph> > long_reads_libs;
    ConvertLongReads(gp.single_long_reads, long_reads_libs);

    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
        io::LibraryType type = cfg::get().ds.reads[i].type();
        if (cfg::get().ds.reads[i].data().mean_insert_size != 0.0 &&
                (type == io::LibraryType::PairedEnd ||
                 type == io::LibraryType::MatePairs)) {

            pe_indexes.push_back(&gp.clustered_indices[i]);
            pe_scaf_indices.push_back(&gp.scaffolding_indices[i]);
            indexes.push_back(i);
        }
    }

#if 0
    if (cfg::get().coverage_based_rr_on == true) {
        std::vector<PathInfo<Graph> > filteredPaths;
        size_t pe_lib_index = GetFirstPELibIndex();
        const io::SequencingLibrary<debruijn_config::DataSetData> &lib =
                cfg::get().ds.reads[pe_lib_index];
        resolve_repeats_by_coverage(gp, (size_t) lib.data().mean_insert_size,
                                    filteredPaths, clustered_indices[0],
                                    quality_labeler);
        PathStorageInfo<Graph> single_storage(
                filteredPaths,
                cfg::get().pe_params.long_reads.coverage_base_rr.filtering,
                cfg::get().pe_params.long_reads.coverage_base_rr.weight_priority,
                cfg::get().pe_params.long_reads.coverage_base_rr.unique_edge_priority);
        long_reads_libs.push_back(single_storage);
    }
#endif

    std::string name = "scaffolds.fasta";
    bool traverse_loops = true;
    if (!(cfg::get().use_scaffolder && cfg::get().pe_params.param_set.scaffolder_options.on)) {
        name = "final_contigs.fasta";
        pe_scaf_indices.clear();
        traverse_loops = false;
    }
    path_extend::ResolveRepeatsPe(
            gp, pe_indexes, pe_scaf_indices, indexes, long_reads_libs,
            cfg::get().output_dir, name, traverse_loops,
            boost::optional<std::string>("final_contigs.fasta"));
}

void RepeatResolution::run(conj_graph_pack &gp) {
    OutputContigs(gp.g, cfg::get().additional_contigs, cfg::get().use_unipaths,
                  cfg::get().simp.tec.plausibility_length);
    OutputContigs(gp.g, cfg::get().output_dir + "before_rr.fasta");

    //What is this?
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

    //todo refactor labeler creation -- and what's that?
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

    if (cfg::get().paired_mode && no_valid_libs && !cfg::get().long_single_mode)
        WARN("Insert size was not estimated for any of the paired libraries, repeat resolution module will not run.");

    if ((no_valid_libs ||
            cfg::get().rm == debruijn_graph::resolving_mode::rm_none) &&
            !cfg::get().long_single_mode) {
        OutputContigs(gp.g, cfg::get().output_dir + "final_contigs.fasta");
    }

    // Repeat resolving begins
    if (cfg::get().rm == debruijn_graph::resolving_mode::rm_path_extend) {
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
