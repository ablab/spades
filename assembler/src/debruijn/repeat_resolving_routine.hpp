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
#include "mismatch_masker.hpp"
#include "contig_output.hpp"

#include "pac_index.hpp"
#include "long_read_storage.hpp"
#include "loop_filter.hpp"
#include "graphio.hpp"
#include "coverage_based_rr.hpp"
#include "pacbio_aligner.hpp"

namespace debruijn_graph {

bool prepare_scaffolding_index(conj_graph_pack& gp,
                               const io::SequencingLibrary<debruijn_config::DataSetData> &lib,
                               PairedIndexT& paired_index,
                               PairedIndexT& clustered_index) {

	double is_var = lib.data().insert_size_deviation;
	size_t delta = size_t(is_var);
	size_t linkage_distance = size_t(
			cfg::get().de.linkage_distance_coeff * is_var);
	GraphDistanceFinder<Graph> dist_finder(gp.g, (size_t) math::round(lib.data().mean_insert_size),
	                                       lib.data().read_length, delta);
	size_t max_distance = size_t(cfg::get().de.max_distance_coeff * is_var);
	boost::function<double(int)> weight_function;

	DEBUG("Retaining insert size distribution for it");
	if (lib.data().insert_size_distribution.size() == 0) {
	    return false;
	}

	WeightDEWrapper wrapper(lib.data().insert_size_distribution, lib.data().mean_insert_size);
	DEBUG("Weight Wrapper Done");
	weight_function = boost::bind(&WeightDEWrapper::CountWeight, wrapper, _1);

	PairInfoWeightFilter<Graph> filter(gp.g, 0.);
	DEBUG("Weight Filter Done");

	const AbstractDistanceEstimator<Graph>& estimator =
			SmoothingDistanceEstimator<Graph>(gp.g, paired_index, dist_finder,
					weight_function, linkage_distance, max_distance,
					cfg::get().ade.threshold, cfg::get().ade.range_coeff,
					cfg::get().ade.delta_coeff, cfg::get().ade.cutoff,
					cfg::get().ade.min_peak_points, cfg::get().ade.inv_density,
					cfg::get().ade.percentage,
					cfg::get().ade.derivative_threshold, true);
	estimate_with_estimator(gp.g, estimator, filter, clustered_index);

	return true;
}

void resolve_repeats_by_coverage(conj_graph_pack& conj_gp, size_t insert_size, std::vector< PathInfo<Graph> >& filteredPaths,
                                 PairedIndexT& clustered_index,
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

	FlankingCoverage<Graph, KmerIndex> index(conj_gp.g, kmer_index, 50);
	EdgeLabelHandler<conj_graph_pack::graph_t> labels_after(conj_gp.g, conj_gp.g);
	auto cov_rr = CoverageBasedResolution<conj_graph_pack> (&conj_gp, cfg::get().cbrr.coverage_threshold_one_list, cfg::get().cbrr.coverage_threshold_match,
			cfg::get().cbrr.coverage_threshold_global, cfg::get().cbrr.tandem_ratio_lower_threshold, cfg::get().cbrr.tandem_ratio_upper_threshold, cfg::get().cbrr.repeat_length_upper_threshold);
	cov_rr.resolve_repeats_by_coverage(index, insert_size, labels_after, quality_labeler, clustered_index, filteredPaths);

	INFO("Repeats are resolved by coverage");
}

size_t get_first_pe_lib_index() {
	for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i)
		if (cfg::get().ds.reads[i].type() == io::LibraryType::PairedEnd && cfg::get().ds.reads[i].data().mean_insert_size != 0.0)
			return i;

    return -1UL;
}

void prepare_all_scaf_libs(conj_graph_pack& conj_gp,
                           vector<PairedIndexT*>& scaff_indexs, vector<size_t>& indexes) {
	vector<PairedIndexT*> cl_scaff_indexs;
	for (size_t i = 0; i < scaff_indexs.size(); ++i) {
		PairedIndexT* pe = new PairedIndexT(conj_gp.g);
		cl_scaff_indexs.push_back(pe);
		INFO("Scaffolding distance estimating started for lib #" << indexes[i]);
		if (!prepare_scaffolding_index(conj_gp, cfg::get().ds.reads[indexes[i]], *scaff_indexs[i], *cl_scaff_indexs[i])) {
		    WARN("Lib #" << indexes[i] << " will not be used for scaffolding");
		}
	}
	scaff_indexs.clear();
	scaff_indexs.insert(scaff_indexs.end(), cl_scaff_indexs.begin(),
			cl_scaff_indexs.end());
}

//todo check usages!!!
void delete_index(vector<PairedIndexT*>& index) {
	for (size_t i = 0; i < index.size(); ++i)
		delete index[i];
}

void pe_resolving(conj_graph_pack& conj_gp, PairedIndicesT& paired_indices,	PairedIndicesT& clustered_indices, const EdgeQuality<Graph, Index>& quality_labeler) {
	vector<PairedIndexT*> pe_indexs;
	vector<PairedIndexT*> pe_scaf_indexs;
	vector<size_t> indexes;

	for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
        const auto& lib = cfg::get().ds.reads[i];
		if (lib.data().mean_insert_size != 0.0 &&
            (lib.type() == io::LibraryType::PairedEnd ||
             lib.type() == io::LibraryType::MatePairs)) {
			pe_indexs.push_back(&clustered_indices[i]);
			pe_scaf_indexs.push_back(&paired_indices[i]);
			indexes.push_back(i);
		}
	}

    //LongReadStorage<Graph> long_read(conj_gp.g);
     //long_read.LoadFromFile("/storage/labnas/students/igorbunova/path-extend2/algorithmic-biology/assembler/pacbio.mpr");

    PathStorage<Graph> long_read(conj_gp.g);
    GapStorage<Graph> gaps(conj_gp.g);

	std::vector< PathInfo<Graph> > filteredPaths;
	if (cfg::get().developer_mode) {
	    OutputContigs(conj_gp.g, cfg::get().output_dir + "before_resolve.fasta");
	}

	if (cfg::get().coverage_based_rr_on == true){
		size_t pe_lib_index = get_first_pe_lib_index();
		const io::SequencingLibrary<debruijn_config::DataSetData> &lib = cfg::get().ds.reads[pe_lib_index];
		resolve_repeats_by_coverage(conj_gp, (size_t) lib.data().mean_insert_size, filteredPaths, clustered_indices[0], quality_labeler);
	}

    //LongReadStorage<Graph> long_read(conj_gp.g);
	if (cfg::get().pacbio_test_on) {
		INFO("creating multiindex with k = " << cfg::get().pb.pacbio_k);
		PacBioAligner pac_aligner(conj_gp, cfg::get().pb.pacbio_k);
		INFO("index created");
		filteredPaths = long_read.GetAllPaths();
		pac_aligner.pacbio_test(long_read, gaps);
	}

	if (cfg::get().use_scaffolder && cfg::get().pe_params.param_set.scaffolder_options.on) {
	    if (cfg::get().pe_params.param_set.scaffolder_options.cluster_info) {
	        prepare_all_scaf_libs(conj_gp, pe_scaf_indexs, indexes);
	    }
        //path_extend::resolve_repeats_pe(conj_gp, pe_indexs, pe_scaf_indexs, indexes, long_read.GetAllPaths(), cfg::get().output_dir, "scaffolds.fasta", true, boost::optional<std::string>("final_contigs.fasta"));
        path_extend::resolve_repeats_pe(conj_gp, pe_indexs, pe_scaf_indexs, indexes, filteredPaths, cfg::get().output_dir, "scaffolds.fasta", true, boost::optional<std::string>("final_contigs.fasta"));
        delete_index(pe_scaf_indexs);
	} else {
		pe_scaf_indexs.clear();
		//path_extend::resolve_repeats_pe(conj_gp, pe_indexs, pe_scaf_indexs, indexes, long_read.GetAllPaths(), cfg::get().output_dir, "final_contigs.fasta", false, boost::none);
		path_extend::resolve_repeats_pe(conj_gp, pe_indexs, pe_scaf_indexs, indexes, filteredPaths, cfg::get().output_dir, "final_contigs.fasta", false, boost::none);
	}
}

void resolve_repeats() {
	Sequence genome = cfg::get().developer_mode ? cfg::get().ds.reference_genome : Sequence();

	conj_graph_pack conj_gp(cfg::get().K, cfg::get().output_dir, genome,
			cfg::get().pos.max_single_gap, cfg::get().pos.careful_labeling,
			!cfg::get().developer_mode);

	PairedIndicesT paired_indices(conj_gp.g, cfg::get().ds.reads.lib_count());
	PairedIndicesT clustered_indices(conj_gp.g,	cfg::get().ds.reads.lib_count());

	if (!cfg::get().developer_mode) {
		conj_gp.edge_pos.Detach();
		paired_indices.Detach();
		clustered_indices.Detach();
		if (!cfg::get().gap_closer_enable && !cfg::get().paired_mode) {
		    //todo ?
//			conj_gp.kmer_mapper.Detach();
		}
	}

	exec_distance_estimation(conj_gp, paired_indices, clustered_indices);

	if (cfg::get().developer_mode && cfg::get().pos.late_threading) {
		FillPos(conj_gp, conj_gp.genome, "10");
		FillPos(conj_gp, !conj_gp.genome, "11");
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

	//todo refactor labeler creation
	total_labeler_graph_struct graph_struct(conj_gp.g, &conj_gp.int_ids,
			&conj_gp.edge_pos);
	total_labeler tot_lab(&graph_struct);
	EdgeQuality<Graph, Index> quality_labeler(conj_gp.g, conj_gp.index,
			conj_gp.kmer_mapper, conj_gp.genome);
	CompositeLabeler<Graph> labeler(tot_lab, quality_labeler);
	detail_info_printer printer(conj_gp, labeler, cfg::get().output_dir);
	printer(ipp_before_repeat_resolution);

	bool no_valid_libs = true;
    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
        if (cfg::get().ds.reads[i].data().mean_insert_size != 0.0) {
            no_valid_libs = false;
            break;
        }
    }

    if (cfg::get().paired_mode && no_valid_libs)
        WARN("Insert size was not estimated for any of the paired libraries, repeat resolution module will not run.");

	if (!cfg::get().paired_mode ||
        no_valid_libs ||
        cfg::get().rm == debruijn_graph::resolving_mode::rm_none) {
		OutputContigs(conj_gp.g, cfg::get().output_dir + "final_contigs.fasta");

        if (cfg::get().pacbio_test_on) {
		    PathStorage<Graph> long_read(conj_gp.g);
		    GapStorage<Graph> gaps(conj_gp.g);
			std::vector< PathInfo<Graph> > filteredPaths;
		    //LongReadStorage<Graph> long_read(conj_gp.g);
			INFO("creating  multiindex with k = " << cfg::get().pb.pacbio_k);
			PacBioAligner pac_aligner(conj_gp, cfg::get().pb.pacbio_k);
			INFO("index created");
			filteredPaths = long_read.GetAllPaths();
			pac_aligner.pacbio_test(long_read, gaps);
		}
		return;
	}

    OutputContigs(conj_gp.g, cfg::get().output_dir + "before_rr.fasta");

	// Repeat resolving begins
	size_t pe_lib_index = get_first_pe_lib_index();
	INFO("STAGE == Resolving Repeats");
    if (cfg::get().ds.reads.lib_count() > 1 || pe_lib_index == -1UL ||
        cfg::get().rm == debruijn_graph::resolving_mode::rm_path_extend) {
		INFO("Path-Extend repeat resolving");
		pe_resolving(conj_gp, paired_indices, clustered_indices, quality_labeler);
	} else if (cfg::get().rm == debruijn_graph::resolving_mode::rm_rectangles) {
		INFO("Ready to run rectangles repeat resolution module");
	} else {
		INFO("Unsupported repeat resolver");
		OutputContigs(conj_gp.g, cfg::get().output_dir + "final_contigs.fasta");
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

