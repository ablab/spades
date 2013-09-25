//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * repeat_resolving_routine.hpp
 *
 *  Created on: 1 Sep 2011
 *      Author: valery
 *      Poor Valery
 */

#pragma once

#include "standard.hpp"

#include "logger/logger.hpp"
#include "distance_estimation_routine.hpp"
#include "io/careful_filtering_reader_wrapper.hpp"
#include "io/is_corrupting_wrapper.hpp"
#include "graph_construction.hpp"
#include "debruijn_stats.hpp"
#include "de/distance_estimation.hpp"
#include "omni/omni_utils.hpp"

#include "path_utils.hpp"
#include "pair_info_improver.hpp"

#include "path_extend/path_extend_launch.hpp"
#include "contig_output.hpp"

#include "pac_index.hpp"
#include "long_read_storage.hpp"
#include "loop_filter.hpp"
#include "graphio.hpp"
#include "coverage_based_rr.hpp"
#include "pacbio_aligner.hpp"
#include "bucket_mapper.hpp"
#include "path_extend/long_read_mapper.hpp"

namespace debruijn_graph {


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
		INFO("Updating index from graph started");

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

/*
 * Return index of first paired-end library or -1 if there is no paired end library
 */
size_t get_first_pe_lib_index() {
    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i)
        if (cfg::get().ds.reads[i].type() == io::LibraryType::PairedEnd && cfg::get().ds.reads[i].data().mean_insert_size != 0.0)
            return i;

    return -1UL;
}

void AddSingleLongReads(vector<PathStorageInfo<Graph> > &long_reads_libs,
                        const vector<PathStorage<Graph>* >& single_long_reads) {
    for (size_t i = 0; i < single_long_reads.size(); ++i) {
        PathStorage<Graph>* storage = single_long_reads[i];
        vector<PathInfo<Graph> > paths = storage->GetAllPaths();
        PathStorageInfo<Graph> single_storage(
                paths,
                cfg::get().pe_params.long_reads.single_reads.filtering,
                cfg::get().pe_params.long_reads.single_reads.weight_priority,
                cfg::get().pe_params.long_reads.single_reads
                        .unique_edge_priority);
        long_reads_libs.push_back(single_storage);
    }
}

void pe_resolving(conj_graph_pack& conj_gp, PairedIndicesT& paired_indexes,
                  PairedIndicesT& clustered_indices,
                  PairedIndicesT& scaffold_indices,
                  const EdgeQuality<Graph, Index>& quality_labeler,
                  vector<PathStorageInfo<Graph> > &long_reads_libs,
                  vector<PathStorage<Graph>* >& single_long_reads) {

    vector<PairedIndexT*> pe_indexes;
    vector<PairedIndexT*> pe_scaf_indices;
    vector<size_t> indexes;
    GapStorage<Graph> gaps(conj_gp.g);
    AddSingleLongReads(long_reads_libs, single_long_reads);

    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
        io::LibraryType type = cfg::get().ds.reads[i].type();
        if (cfg::get().ds.reads[i].data().mean_insert_size != 0.0 &&
                (type == io::LibraryType::PairedEnd
                || type == io::LibraryType::MatePairs)) {

            pe_indexes.push_back(&clustered_indices[i]);
            pe_scaf_indices.push_back(&scaffold_indices[i]);
            indexes.push_back(i);
        }
    }

    if (cfg::get().coverage_based_rr_on == true) {
        std::vector<PathInfo<Graph> > filteredPaths;
        size_t pe_lib_index = get_first_pe_lib_index();
        const io::SequencingLibrary<debruijn_config::DataSetData> &lib =
                cfg::get().ds.reads[pe_lib_index];
        resolve_repeats_by_coverage(conj_gp, (size_t) lib.data().mean_insert_size,
                                    filteredPaths, clustered_indices[0],
                                    quality_labeler);
        PathStorageInfo<Graph> single_storage(
                filteredPaths,
                cfg::get().pe_params.long_reads.coverage_base_rr.filtering,
                cfg::get().pe_params.long_reads.coverage_base_rr.weight_priority,
                cfg::get().pe_params.long_reads.coverage_base_rr.unique_edge_priority);
        long_reads_libs.push_back(single_storage);
    }

    std::string name = "scaffolds.fasta";
    bool traverse_loops = true;
    if (!(cfg::get().use_scaffolder
            && cfg::get().pe_params.param_set.scaffolder_options.on)) {
        name = "final_contigs.fasta";
        pe_scaf_indices.clear();
        traverse_loops = false;
    }
    path_extend::ResolveRepeatsPe(
            conj_gp, pe_indexes, pe_scaf_indices, indexes, long_reads_libs,
            cfg::get().output_dir, name, traverse_loops,
            boost::optional<std::string>("final_contigs.fasta"));
}

void resolve_repeats() {
	Sequence genome = cfg::get().developer_mode ? cfg::get().ds.reference_genome : Sequence();

	conj_graph_pack conj_gp(cfg::get().K, cfg::get().output_dir, genome,
			cfg::get().pos.max_single_gap, cfg::get().pos.careful_labeling,
			!cfg::get().developer_mode);

	PairedIndicesT paired_indices(conj_gp.g, cfg::get().ds.reads.lib_count());
	PairedIndicesT clustered_indices(conj_gp.g,	cfg::get().ds.reads.lib_count());
    PairedIndicesT scaffold_indices(conj_gp.g, cfg::get().ds.reads.lib_count());
    vector<PathStorageInfo<Graph> > long_reads_libs;
    vector<PathStorage<Graph>* > single_long_reads;

    for (size_t i = 0; i < cfg::get().ds.count_single_libs; ++i) {
        single_long_reads.push_back(new PathStorage<Graph>(conj_gp.g));
    }
	if (!cfg::get().developer_mode) {
		conj_gp.edge_pos.Detach();
		paired_indices.Detach();
		clustered_indices.Detach();
		if (!cfg::get().gap_closer_enable && !cfg::get().paired_mode) {
		    //todo ?
//			conj_gp.kmer_mapper.Detach();
		}
	 }
     PathStorage<Graph> pacbio_read(conj_gp.g);

	 exec_distance_estimation(conj_gp, paired_indices, clustered_indices, scaffold_indices, pacbio_read, single_long_reads);

	 if (cfg::get().entry_point <= ws_pacbio_aligning && cfg::get().gap_closer_enable){
	     INFO(" need to align pb");
	     if (cfg::get().pacbio_test_on) {
	         INFO("creating  multiindex with k = " << cfg::get().pb.pacbio_k);
	         PacBioAligner pac_aligner(conj_gp, paired_indices, clustered_indices, scaffold_indices, single_long_reads, cfg::get().pb.pacbio_k);
	         INFO("index created");
	         GapStorage<Graph> gaps(conj_gp.g);
	         pac_aligner.pacbio_test(pacbio_read, gaps);
	     }
         for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
             io::LibraryType type = cfg::get().ds.reads[i].type();
             if (type == io::LibraryType::PacBioReads) {
                 //TODO: need to read reads from stream instead of file and delete pacbio_on + pacbio reads from config
                 PathStorage<Graph> pacbio_read1(conj_gp.g);
                 INFO("creating  multiindex with k = " << cfg::get().pb.pacbio_k);
                 PacBioAligner pac_aligner(conj_gp, paired_indices, clustered_indices, scaffold_indices, single_long_reads, cfg::get().pb.pacbio_k);
                 INFO("index created");
                 GapStorage<Graph> gaps(conj_gp.g);
                 pac_aligner.pacbio_test(pacbio_read1, gaps);
                 vector<PathInfo<Graph> > pacbio_paths = pacbio_read.GetAllPaths();
                 PathStorageInfo<Graph> pacbio_storage(
                         pacbio_paths,
                         cfg::get().pe_params.long_reads.pacbio_reads.filtering,
                         cfg::get().pe_params.long_reads.pacbio_reads.weight_priority,
                         cfg::get().pe_params.long_reads.pacbio_reads.unique_edge_priority);
                 long_reads_libs.push_back(pacbio_storage);
             }
        }
	}
	if (cfg::get().pacbio_test_on && cfg::get().gap_closer_enable) {
	    INFO("getting paths");
        vector<PathInfo<Graph> > pacbio_paths = pacbio_read.GetAllPaths();
        PathStorageInfo<Graph> pacbio_storage(
        pacbio_paths,
        cfg::get().pe_params.long_reads.pacbio_reads.filtering,
        cfg::get().pe_params.long_reads.pacbio_reads.weight_priority,
        cfg::get().pe_params.long_reads.pacbio_reads.unique_edge_priority);
        INFO("storage created");
        long_reads_libs.push_back(pacbio_storage);
	}
	if (cfg::get().developer_mode && cfg::get().pos.late_threading) {
	    INFO("threading");
		FillPos(conj_gp, conj_gp.genome, "10");
		FillPos(conj_gp, !conj_gp.genome, "11");
		INFO("and");
		if (!cfg::get().pos.contigs_for_threading.empty()
				&& FileExists(cfg::get().pos.contigs_for_threading)) {
			FillPosWithRC(conj_gp, cfg::get().pos.contigs_for_threading,
					"thr_");
		}

		if (!cfg::get().pos.contigs_to_analyze.empty()
				&& FileExists(cfg::get().pos.contigs_to_analyze)) {
			FillPosWithRC(conj_gp, cfg::get().pos.contigs_to_analyze, "anlz_");
		}
	}


//	RunTopologyTipClipper(conj_gp.g, 300, 2000, 1000);
	INFO("threaded");
	//todo refactor labeler creation
	total_labeler_graph_struct graph_struct(conj_gp.g, &conj_gp.int_ids,
			&conj_gp.edge_pos);
	total_labeler tot_lab(&graph_struct);
	EdgeQuality<Graph, Index> quality_labeler(conj_gp.g, conj_gp.index,
			conj_gp.kmer_mapper, conj_gp.genome);
	CompositeLabeler<Graph> labeler(tot_lab, quality_labeler);
	INFO("lavelers created");
	detail_info_printer printer(conj_gp, labeler, cfg::get().output_dir);
	INFO("printing");
	printer(ipp_before_repeat_resolution);

	bool no_valid_libs = true;
    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
        if (cfg::get().ds.reads[i].data().mean_insert_size != 0.0) {
            no_valid_libs = false;
            break;
        }
    }

    if (cfg::get().paired_mode && no_valid_libs && !cfg::get().long_single_mode && !cfg::get().pacbio_test_on) {
        WARN("Insert size was not estimated for any of the paired libraries, repeat resolution module will not run.");
    }

	if ((!cfg::get().paired_mode
	        || no_valid_libs
			|| cfg::get().rm == debruijn_graph::resolving_mode::rm_none) && !cfg::get().long_single_mode && !cfg::get().pacbio_test_on) {
		OutputContigs(conj_gp.g, cfg::get().output_dir + "final_contigs.fasta");
		return;
	}

    OutputContigs(conj_gp.g, cfg::get().output_dir + "before_rr.fasta");

	//Repeat resolving begins
	size_t pe_lib_index = get_first_pe_lib_index();
	INFO("STAGE == Resolving Repeats");
	if (cfg::get().long_single_mode || cfg::get().ds.reads.lib_count() > 1 || pe_lib_index == -1UL
			|| cfg::get().rm
					== debruijn_graph::resolving_mode::rm_path_extend) {
		INFO("Path-Extend repeat resolving");
		pe_resolving(conj_gp, paired_indices, clustered_indices,  scaffold_indices, quality_labeler, long_reads_libs, single_long_reads);
	}
	else if (cfg::get().rm == debruijn_graph::resolving_mode::rm_rectangles) {
		INFO("Ready to run rectangles repeat resolution module");
	} else {
		INFO("Unsupported repeat resolver");
		OutputContigs(conj_gp.g, cfg::get().output_dir + "final_contigs.fasta");
	}
	for (size_t i = 0; i < cfg::get().ds.count_single_libs; ++i) {
        delete single_long_reads[i];
    }
}

void exec_repeat_resolving() {
	if (cfg::get().entry_point <= ws_repeats_resolving) {
		resolve_repeats();
		//todo why nothing to save???
		// nothsng to save yet
	} else {
		INFO("Loading Repeat Resolving");
		INFO("Nothing to load");
		// nothing to load
	}
}

} // debruijn_graph

