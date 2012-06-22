//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * repeat_resolving_routine.hpp
 *
 *  Created on: 1 Sep 2011
 *      Author: valery
 */

#pragma once

#include "standard.hpp"

#include "logger/logger.hpp"
#include "repeat_resolving.hpp"
#include "distance_estimation_routine.hpp"
//#include "path_set_graph_constructor.hpp"
#include "io/careful_filtering_reader_wrapper.hpp"
#include "io/is_corrupting_wrapper.hpp"
#include "resolved_pair_info.hpp"
#include "graph_construction.hpp"
#include "debruijn_stats.hpp"
#include "omni/distance_estimation.hpp"
#include "omni/omni_utils.hpp"
#include "long_contigs/lc_launch.hpp"
#include "internal_aligner.hpp"
#include "omni/loop_killer.hpp"
#include "path_utils.hpp"

typedef io::CarefulFilteringReaderWrapper<io::SingleRead> CarefulFilteringStream;

namespace debruijn_graph {

void resolve_repeats(PairedReadStream& stream, const Sequence& genome);
} // debruijn_graph

// TODO move impl to *.cpp

namespace debruijn_graph {

void save_distance_filling(conj_graph_pack& gp, paired_info_index& paired_index,
		paired_info_index& clustered_index) {
	if (cfg::get().make_saves) {
		fs::path p = fs::path(cfg::get().output_saves) / "distance_filling";
		PrintAll(p.string(), gp, paired_index, clustered_index);
		write_estimated_params(p.string());
	}
}
void FillContigNumbers(map<ConjugateDeBruijnGraph::EdgeId, int>& contigNumbers
		, ConjugateDeBruijnGraph& cur_graph) {
	int cur_num = 0;
	set<EdgeId> edges;
	for (auto iter = cur_graph.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
		if (edges.find(*iter) == edges.end()) {
			contigNumbers[*iter] = cur_num;
			cur_num++;
			edges.insert(cur_graph.conjugate(*iter));
		}
	}
}

void FillContigNumbers(
		map<NonconjugateDeBruijnGraph::EdgeId, int>& contigNumbers
		, NonconjugateDeBruijnGraph& cur_graph) {
	int cur_num = 0;
	for (auto iter = cur_graph.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
		contigNumbers[*iter] = cur_num;
		cur_num++;
	}
}

int ContigNumber(map<NonconjugateDeBruijnGraph::EdgeId, int>& contigNumbers
		, NonconjugateDeBruijnGraph::EdgeId eid
		, NonconjugateDeBruijnGraph& cur_graph) {
	if (contigNumbers.find(eid) != contigNumbers.end())
		return (contigNumbers[eid]);
	else {
		WARN("Deleted edge");
		return -1;
	}
}

int ContigNumber(map<ConjugateDeBruijnGraph::EdgeId, int>& contigNumbers
		, ConjugateDeBruijnGraph::EdgeId eid
		, ConjugateDeBruijnGraph& cur_graph) {
	if (contigNumbers.find(eid) != contigNumbers.end())
		return (contigNumbers[eid]);
	else if (contigNumbers.find(cur_graph.conjugate(eid))
			!= contigNumbers.end())
		return (contigNumbers[cur_graph.conjugate(eid)]);
	else {
		WARN("Deleted edge");
		return -1;
	}
}

template<size_t k, class graph_pack>
void SelectReadsForConsensusBefore(graph_pack& etalon_gp,
		typename graph_pack::graph_t& cur_graph,
		EdgeLabelHandler<typename graph_pack::graph_t>& LabelsAfter,
		const EdgeIndex<k + 1, typename graph_pack::graph_t>& index,
		vector<ReadStream*> reads, string& consensus_output_dir) {
	INFO("ReadMapping started");
	map<typename graph_pack::graph_t::EdgeId, int> contigNumbers;
	int cur_num = 0;
	FillContigNumbers(contigNumbers, cur_graph);
	for (auto iter = etalon_gp.g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
		DEBUG(
				"Edge number:" << etalon_gp.int_ids.ReturnIntId(*iter)
						<< " is contained in contigs");
		auto images_it = LabelsAfter.edge_inclusions.find(*iter);
		if (images_it != LabelsAfter.edge_inclusions.end()) {
			for (auto it = images_it->second.begin();
					it != images_it->second.end(); ++it) {
				DEBUG(ContigNumber(contigNumbers, *it, cur_graph) << ", ");
			}
		}
	}
	cur_num = contigNumbers.size();
	DEBUG(cur_num << " contigs");
	for (int i = 1; i < 3; i++) {
		int read_num = 0;
		osequencestream* mapped_reads[5000];
		for (int j = 0; j < cur_num; j++) {
			string output_filename = consensus_output_dir + ToString(j)
					+ "_reads" + ToString(i) + ".fa";
			osequencestream* tmp = new osequencestream(output_filename);
//          mapped_reads.push_back(tmp);
			mapped_reads[j] = tmp;
		}
		SingleReadMapper<k, typename graph_pack::graph_t> rm(etalon_gp.g,
				index);
		DEBUG("mapping reads from pair " << i);
		while (!reads[i - 1]->eof()) {
			io::SingleRead cur_read;

			(*reads[i - 1]) >> cur_read;
			vector<typename graph_pack::graph_t::EdgeId> res =
					rm.GetContainingEdges(cur_read);
			read_num++;
			TRACE(
					read_num << " mapped to" << res.size() << " contigs :, read"
							<< cur_read.sequence());
//          map_quantity += res.size();
			for (size_t ii = 0; ii < res.size(); ii++) {
				TRACE("counting number " << contigNumbers[res[ii]]);
				if (ContigNumber(contigNumbers, res[ii], cur_graph) != -1)
					(*mapped_reads[ContigNumber(contigNumbers, res[ii],
							cur_graph)]) << cur_read.sequence();
				else
					WARN(
							"No edges containing"
									<< etalon_gp.int_ids.ReturnIntId(res[ii]));
			}
		}
	}
}

template<size_t k, class graph_pack>
void SelectReadsForConsensus(graph_pack& etalon_gp,
		typename graph_pack::graph_t& cur_graph,
		EdgeLabelHandler<typename graph_pack::graph_t>& LabelsAfter,
		const EdgeIndex<k + 1, typename graph_pack::graph_t>& index
		,vector<ReadStream *>& reads , string& consensus_output_dir) {
	INFO("ReadMapping started");
	map<typename graph_pack::graph_t::EdgeId, int> contigNumbers;
	int cur_num = 0;
	FillContigNumbers(contigNumbers, cur_graph);
	for (auto iter = etalon_gp.g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
		DEBUG(
				"Edge number:" << etalon_gp.int_ids.ReturnIntId(*iter)
						<< " is contained in contigs");
		set<typename graph_pack::graph_t::EdgeId> images =
				LabelsAfter.edge_inclusions[*iter];
		for (auto it = images.begin(); it != images.end(); ++it) {
			DEBUG(ContigNumber(contigNumbers, *it, cur_graph) << ", ");
		}
	}
	cur_num = contigNumbers.size();
	DEBUG(cur_num << " contigs");
	for (int i = 1; i < 3; i++) {
		int read_num = 0;
		osequencestream* mapped_reads[5000];
		for (int j = 0; j < cur_num; j++) {
			string output_filename = consensus_output_dir + ToString(j)
					+ "_reads" + ToString(i) + ".fa";
			osequencestream* tmp = new osequencestream(output_filename);
//          mapped_reads.push_back(tmp);
			mapped_reads[j] = tmp;
		}
		SingleReadMapper<k, typename graph_pack::graph_t> rm(etalon_gp.g,
				index);
		DEBUG("mapping reads from pair " << i);
		while (!reads[i - 1]->eof()) {
			io::SingleRead cur_read;

			(*reads[i - 1]) >> cur_read;
			vector<typename graph_pack::graph_t::EdgeId> res =
					rm.GetContainingEdges(cur_read);
			read_num++;
			TRACE(
					read_num << " mapped to" << res.size() << " contigs :, read"
							<< cur_read.sequence());
//          map_quantity += res.size();
			for (size_t ii = 0; ii < res.size(); ii++) {
				TRACE("counting number " << contigNumbers[res[ii]]);
				set<typename graph_pack::graph_t::EdgeId> images =
						LabelsAfter.edge_inclusions[res[ii]];
				for (auto iter = images.begin(); iter != images.end(); ++iter)
					if (ContigNumber(contigNumbers, *iter, cur_graph) != -1)
						(*mapped_reads[ContigNumber(contigNumbers, *iter,
								cur_graph)]) << cur_read.sequence();
					else
						WARN(
								"No edges containing"
										<< etalon_gp.int_ids.ReturnIntId(
												res[ii]));
			}
		}
	}
}

void SAM_after_resolve(conj_graph_pack& conj_gp, conj_graph_pack& resolved_gp,
		EdgeLabelHandler<conj_graph_pack::graph_t> &labels_after) {
	if (cfg::get().SAM_writer_enable && cfg::get().sw.align_after_RR) {
		int offset = determine_offset(
				input_file(cfg::get().ds.paired_reads[0][0]));
		io::OffsetType offset_type;
		if (offset == 33) {
			INFO("Using offset +33");
			offset_type = io::PhredOffset;
		} else if (offset == 64) {
			INFO("Using offset +64");
			offset_type = io::SolexaOffset;
		} else {
			WARN("Unable to define offset type, assume +33");
			offset_type = io::PhredOffset;
		}

		if (cfg::get().sw.align_original_reads) {
//			if (cfg::get().sw.original_first && cfg::get().sw.original_second)
			{
				auto paired_reads = paired_easy_reader(false, 0, false, false,
						false, offset_type);
				auto original_paired_reads = paired_easy_reader(false, 0, false,
						false, true, offset_type);
//				io::PairedEasyReader original_paired_reads(
//								make_pair(input_file(*cfg::get().sw.original_first),
//										input_file(*cfg::get().sw.original_second)),
//								false,
//								0);
				typedef NewExtendedSequenceMapper<K + 1, Graph> SequenceMapper;
				SequenceMapper mapper(conj_gp.g, conj_gp.index,
						conj_gp.kmer_mapper);

				bool print_quality = (
						cfg::get().sw.print_quality ?
								*cfg::get().sw.print_quality : false);
				OriginalReadsResolvedInternalAligner<ConjugateDeBruijnGraph,
						SequenceMapper> Aligner(resolved_gp.g, conj_gp.g,
						mapper, labels_after, cfg::get().sw.adjust_align,
						cfg::get().sw.output_map_format,
						cfg::get().sw.output_broken_pairs, print_quality);
				Aligner.AlignPairedReads(*original_paired_reads, *paired_reads,
						cfg::get().output_dir + "align_after_RR.sam");
			}
		} else {
			auto paired_reads = paired_easy_reader(false, 0, false, false,
					false, offset_type);
			auto single_reads = single_easy_reader(false, false, false,
					offset_type);

			typedef NewExtendedSequenceMapper<K + 1, Graph> SequenceMapper;
			SequenceMapper mapper(conj_gp.g, conj_gp.index,
					conj_gp.kmer_mapper);

			bool print_quality = (
					cfg::get().sw.print_quality ?
							*cfg::get().sw.print_quality : false);
			ResolvedInternalAligner<ConjugateDeBruijnGraph, SequenceMapper> Aligner(
					resolved_gp.g, conj_gp.g, mapper, labels_after,
					cfg::get().sw.adjust_align, cfg::get().sw.output_map_format,
					cfg::get().sw.output_broken_pairs, print_quality);
			if (cfg::get().sw.align_only_paired)
				Aligner.AlignPairedReads(*paired_reads,
						cfg::get().output_dir + "align_after_RR.sam");
			else
				Aligner.AlignReads(*paired_reads, *single_reads,
						cfg::get().output_dir + "align_after_RR.sam");

		}
	}
}

template<class graph_pack>
void CleanIsolated(graph_pack& gp) {
	for (auto iter = gp.g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
		typename graph_pack::graph_t::VertexId start, end;
		start = gp.g.EdgeStart(*iter);
		end = gp.g.EdgeEnd(*iter);
		TRACE(
				gp.g.CheckUniqueOutgoingEdge(start) << " "
						<< gp.g.IsDeadStart(start) << " "
						<< gp.g.CheckUniqueIncomingEdge(end) << " "
						<< gp.g.IsDeadEnd(end));
		if (gp.g.CheckUniqueOutgoingEdge(start) && gp.g.IsDeadStart(start)
				&& gp.g.CheckUniqueIncomingEdge(end) && gp.g.IsDeadEnd(end))
			gp.g.DeleteEdge(*iter);
	}
}

string GeneratePostfix() {
	string s = "_";
	if (cfg::get().rr.symmetric_resolve)
		s += "sym_";
	else
		s += "nonsym_";

	if (cfg::get().advanced_estimator_mode)
		s += "advanced_est_";
	else
		s += "usual_est_";

	s += "k";
	s += ToString(K);
	if (cfg::get().path_set_graph) {
		s += "_path_set";
	}
	if (cfg::get().rr.mode == 2) {
		s += "_mode2";
	} else {
		s += "_nv";
		s += ToString(cfg::get().rr.near_vertex);
	}
	s += ".fasta";
	return s;
}

template<class graph_pack>
void ProduceResolvedPairedInfo(
		graph_pack& origin_gp,
		PairedInfoIndex<typename graph_pack::graph_t>& clustered_index,
		graph_pack& resolved_gp,
		EdgeLabelHandler<typename graph_pack::graph_t>& labels_after,
		PairedInfoIndex<typename graph_pack::graph_t>& resolved_graph_paired_info) {
	INFO("Generating paired info for resolved graph");
	ResolvedGraphPairInfoCounter<typename graph_pack::graph_t> resolved_graph_paired_info_counter(
			origin_gp.g, clustered_index, resolved_gp.g, labels_after);
	resolved_graph_paired_info_counter.FillResolvedGraphPairedInfo(
			resolved_graph_paired_info);
	DEBUG("Generating paired info for resolved graph complete");
}

template<class graph_pack>
void SaveResolvedPairedInfo(
		graph_pack& resolved_gp,
		PairedInfoIndex<typename graph_pack::graph_t> resolved_graph_paired_info,
		const string& graph_name, const string& subfolder) {
	if (cfg::get().make_saves) {
		std::string rr_filename;
		if (subfolder.size()) {
			INFO("Saving graph and paired info to subfolder " << subfolder);
			rr_filename = (cfg::get().output_dir + subfolder) + graph_name;
		} else {
			INFO("Saving graph and paired info");
			fs::path p = fs::path(cfg::get().output_saves) / (graph_name);
			rr_filename = p.string();
		}
		PrintWithClusteredIndex(rr_filename, resolved_gp,
				resolved_graph_paired_info);
		DEBUG("Saved");
	}
}

template<class graph_pack>
void process_resolve_repeats(graph_pack& origin_gp,
		PairedInfoIndex<typename graph_pack::graph_t>& clustered_index,
		graph_pack& resolved_gp, const string& graph_name,
		EdgeLabelHandler<typename graph_pack::graph_t>& labels_after,
		const string& subfolder = "", bool output_contigs = true) {

//	EdgeLabelHandler<typename graph_pack::graph_t> labels_after(resolved_gp.g,
//			origin_gp.g);
//	ProduceLongEdgesStat( origin_gp,  clustered_index);
	if (cfg::get().compute_paths_number)
		GenerateMatePairStats(origin_gp, clustered_index);
	DEBUG("New index size: " << clustered_index.size());
	// todo: make printGraph const to its arguments

	// todo: possibly we don't need it
//    if (cfg::get().rectangle_mode)
//        RectangleResolve(clustered_index, origin_gp.g, cfg::get().output_root + "tmp/", cfg::get().output_dir);
	string postfix = GeneratePostfix();
	typedef TotalLabelerGraphStruct<typename graph_pack::graph_t> total_labeler_gs;
	typedef TotalLabeler<typename graph_pack::graph_t> total_labeler;
	total_labeler_gs graph_struct_before(origin_gp.g, &origin_gp.int_ids,
			&origin_gp.edge_pos, NULL);
	total_labeler tot_labeler_before(&graph_struct_before);
	if (cfg::get().path_set_graph) {
		VERIFY(false);
//		INFO("testing path-set graphs");
//		PathSetGraphConstructor<graph_pack> path_set_constructor(origin_gp,
//				clustered_index, resolved_gp);
//		path_set_constructor.Construct();
//		typedef typename graph_pack::graph_t::EdgeId EdgeId;
//		map<EdgeId, EdgeId> edge_labels;
//		map<EdgeId, EdgeId> rmap(path_set_constructor.GetEdgeLabels());
//		rmap.Copy(edge_labels);
//		labels_after.FillLabels(edge_labels);
//		INFO("testing ended");
	} else {

//    CleanIsolated(origin_gp);
		ResolveRepeats(
				origin_gp.g,
				origin_gp.int_ids,
				clustered_index,
				origin_gp.edge_pos,
				resolved_gp.g,
				resolved_gp.int_ids,
				resolved_gp.edge_pos,
				cfg::get().output_dir + subfolder + "resolve_" + graph_name
						+ "/", labels_after);

		//Generating paired info for resolved graph
//		PairedInfoIndex<typename graph_pack::graph_t> resolved_graph_paired_info(resolved_gp.g);
//		ProduceResolvedPairedInfo(origin_gp, clustered_index, resolved_gp, labels_after, resolved_graph_paired_info);
//		SaveResolvedPairedInfo(resolved_gp, resolved_graph_paired_info, graph_name + "_resolved", subfolder);
		//Paired info for resolved graph generated

	}

	if (cfg::get().output_nonfinal_contigs && output_contigs) {
		OutputContigs(resolved_gp.g,
				cfg::get().output_dir + "after_rr_before_simplify" + postfix);
//		OutputContigs(origin_gp.g,
//				cfg::get().output_dir + "before_resolve" + postfix);
	}
	INFO("Running total labeler");

	total_labeler_gs graph_struct_after(resolved_gp.g, &resolved_gp.int_ids,
			&resolved_gp.edge_pos, &labels_after);
	total_labeler tot_labeler_after(&graph_struct_after, &graph_struct_before);
	if (cfg::get().output_pictures) {
		omnigraph::WriteSimple(
				resolved_gp.g,
				tot_labeler_after,
				cfg::get().output_dir + subfolder + graph_name
						+ "_3_resolved.dot", "no_repeat_graph");
	}

	DEBUG("Total labeler finished");

	//Generating paired info for resolved graph
	{
		PairedInfoIndex<typename graph_pack::graph_t> resolved_cleared_graph_paired_info_before(
				resolved_gp.g);
		if (cfg::get().path_set_graph == false) {

			ProduceResolvedPairedInfo(origin_gp, clustered_index, resolved_gp,
					labels_after, resolved_cleared_graph_paired_info_before);
		}
	}

	INFO("SUBSTAGE == Clearing resolved graph");

	omnigraph::Compressor<typename graph_pack::graph_t> compressor(
			resolved_gp.g);
	compressor.CompressAllVertices();

//	    omnigraph::StrGraphLabeler<typename graph_pack::graph_t> str_labeler(resolved_gp.g);
//	omnigraph::WriteSimple(resolved_gp.g, str_labeler,
//				cfg::get().output_dir + subfolder + graph_name + "_resolved.dot",
//				"no_repeat_graph");

	EdgeRemover<typename graph_pack::graph_t> edge_remover(resolved_gp.g,
			false);
	size_t iters = 3; // TODO Constant 3? Shouldn't it be taken from config?
	for (size_t i = 0; i < iters; ++i) {
		INFO(
				"Tip clipping iteration " << i << " (0-indexed) out of "
						<< iters << ":");

		ClipTipsForResolver(resolved_gp.g);

		if (cfg::get().path_set_graph == false) {
			PairedInfoIndex<typename graph_pack::graph_t> resolved_cleared_graph_paired_info_before(
					resolved_gp.g);

			ProduceResolvedPairedInfo(origin_gp, clustered_index, resolved_gp,
					labels_after, resolved_cleared_graph_paired_info_before);
		}

//		INFO("Erroneous remove "<<i);
//        BulgeRemoveWrap      (resolved_gp.g);
//		FinalRemoveErroneousEdges(resolved_gp.g, edge_remover);

		if (cfg::get().path_set_graph == false) {
			PairedInfoIndex<typename graph_pack::graph_t> resolved_cleared_graph_paired_info_before(
					resolved_gp.g);

			ProduceResolvedPairedInfo(origin_gp, clustered_index, resolved_gp,
					labels_after, resolved_cleared_graph_paired_info_before);
		}

//        RemoveRelativelyLowCoverageEdges(resolved_gp.g);
//		omnigraph::WriteSimple(resolved_gp.g, tot_labeler_after,
//
//		cfg::get().output_dir + subfolder + ToString(i) + "b_4_cleared.dot",
//				"no_repeat_graph");
	}

//	SimpleLoopKiller<typename graph_pack::graph_t> lk(resolved_gp.g, 3000, 6);
//	lk.KillAllLoops();

	DEBUG("Clearing resolved graph complete");

	//Generating paired info for resolved graph
	if (cfg::get().path_set_graph == false) {
		PairedInfoIndex<typename graph_pack::graph_t> resolved_cleared_graph_paired_info(
				resolved_gp.g);

		ProduceResolvedPairedInfo(origin_gp, clustered_index, resolved_gp,
				labels_after, resolved_cleared_graph_paired_info);
		SaveResolvedPairedInfo(resolved_gp, resolved_cleared_graph_paired_info,
				graph_name + "_resolved_cleared", subfolder);
	}
	//Paired info for resolved graph generated

	DEBUG("Output Contigs");

	if (output_contigs) {
//		if (cfg::get().output_nonfinal_contigs) {
//		    OutputContigs(resolved_gp.g, cfg::get().output_dir + "resolved_and_cleared" + postfix);
//		}
		OutputContigs(resolved_gp.g,
				cfg::get().output_dir + "final_contigs.fasta");
		cfg::get_writable().final_contigs_file = cfg::get().output_dir
				+ "final_contigs.fasta";
	}

	if (cfg::get().output_pictures) {
		omnigraph::WriteSimple(
				resolved_gp.g,
				tot_labeler_after,
				cfg::get().output_dir + subfolder + graph_name
						+ "_4_cleared.dot", "no_repeat_graph");
	}

	if (output_contigs) {
		if (cfg::get().need_consensus) {
			string consensus_folder = cfg::get().output_dir
					+ "consensus_after_resolve/";
			OutputSingleFileContigs(resolved_gp.g, consensus_folder);

//			SelectReadsForConsensus<K,  graph_pack>(origin_gp, resolved_gp.g, labels_after, origin_gp.index, reads, consensus_folder);
			consensus_folder = cfg::get().output_dir
					+ "consensus_before_resolve/";
			OutputSingleFileContigs(origin_gp.g, consensus_folder);
			VERIFY_MSG(
					false,
					"Due to compilation reasons, need_consensus is deprecated.");
//			vector<ReadStream*> streams = single_streams(true, true);
//			SelectReadsForConsensusBefore<K, graph_pack>(origin_gp, origin_gp.g,
//					labels_after, origin_gp.index, streams, consensus_folder);

		}

		/*		one_many_contigs_enlarger<typename graph_pack::graph_t> N50enlarger(
		 resolved_gp.g, cfg::get().ds.IS);
		 N50enlarger.Loops_resolve();

		 omnigraph::WriteSimple(
		 resolved_gp.g,
		 tot_labeler_after,
		 cfg::get().output_dir + subfolder + graph_name
		 + "_5_unlooped.dot", "no_repeat_graph");

		 OutputContigs(resolved_gp.g,
		 cfg::get().output_dir + "unlooped" + postfix);
		 */
//		N50enlarger.one_many_resolve_with_vertex_split();
//
//		omnigraph::WriteSimple(
//				resolved_gp.g,
//				tot_labeler_after,
//				cfg::get().output_dir + subfolder + graph_name
//						+ "_6_finished.dot", "no_repeat_graph");
//
//		OutputContigs(resolved_gp.g,
//				cfg::get().output_dir + "contigs_final.fasta");
	}
}


template<class graph_pack>
bool TryToAddPairInfo(const graph_pack &origin_gp,
		PairedInfoIndex<typename graph_pack::graph_t>& clustered_index,
		const PairInfo<typename graph_pack::graph_t::EdgeId>& p_i) {
		return TryToAddPairInfo(origin_gp, clustered_index, p_i.first, p_i.second, p_i.d, p_i.weight, p_i.variance);
}



template<class graph_pack>
bool TryToAddPairInfo(const graph_pack &origin_gp,
		PairedInfoIndex<typename graph_pack::graph_t>& clustered_index,
		typename graph_pack::graph_t::EdgeId first,
		typename graph_pack::graph_t::EdgeId second, double tmpd, double w, double var = 0) {
	auto pi_vector = clustered_index.GetEdgePairInfo(first, second);
	bool already_exist = false;
	DEBUG(
			"TRY TO adding paired info between edges "
					<< origin_gp.int_ids.ReturnIntId(first) << " "
					<< origin_gp.int_ids.ReturnIntId(second) << " dist " << tmpd
					<< " Size of possible conflicts " << pi_vector.size());

	for (auto pi_iter = pi_vector.begin(); pi_iter != pi_vector.end();
			++pi_iter) {
		DEBUG(
				"dist " << tmpd << " versus " << pi_iter->d << " var "
						<< pi_iter->variance);
		if (abs(pi_iter->d - tmpd) < 0.01 + pi_iter->variance + var) {
			already_exist = true;
		}
	}
	if (!already_exist) {
//		added_info++;
		INFO(
				"adding paired info between edges "
						<< origin_gp.int_ids.ReturnIntId(first) << " "
						<< origin_gp.int_ids.ReturnIntId(second) << " dist "
						<< tmpd);
		clustered_index.AddPairInfo(
				PairInfo<typename graph_pack::graph_t::EdgeId>(first, second,
						tmpd, w, 0));
		INFO(
				"adding paired info between edges "
						<< origin_gp.int_ids.ReturnIntId(
								origin_gp.g.conjugate(second))
						<< " "
						<< origin_gp.int_ids.ReturnIntId(
								origin_gp.g.conjugate(first))
						<< " dist "
						<< tmpd - origin_gp.g.length(first)
								+ origin_gp.g.length(second));
		clustered_index.AddPairInfo(
				PairInfo<typename graph_pack::graph_t::EdgeId>(
						origin_gp.g.conjugate(second),
						origin_gp.g.conjugate(first),
						tmpd - origin_gp.g.length(first)
								+ origin_gp.g.length(second), w, 0));
		return true;
	}
	return false;
}



template<class graph_pack>
int DeleteIfExist(const graph_pack &origin_gp,
		PairedInfoIndex<typename graph_pack::graph_t>& clustered_index,
		const PairInfo<typename graph_pack::graph_t::EdgeId>& p_info) {
	int cnt = 0;
	auto pi_vector = clustered_index.GetEdgePairInfo(p_info.first, p_info.second);
	for (auto pi_iter = pi_vector.begin(); pi_iter != pi_vector.end();
			++pi_iter) {
		if (abs(pi_iter->d - p_info.d) < 0.01) {
			clustered_index.RemovePairInfo(*pi_iter);
			clustered_index.RemovePairInfo(BackwardInfo(*pi_iter));
			cnt += 2;
			INFO("Removed pi "<< origin_gp.g.int_id(pi_iter->first) << " "<< origin_gp.g.int_id(pi_iter->second)<<" d "<<pi_iter->d<<" var "<<pi_iter->variance);
		}
	}
	INFO("cnt += "<<cnt);
	return cnt;
}



template<class graph_pack>
int DeleteConjugatePairInfo(const graph_pack &origin_gp,
		PairedInfoIndex<typename graph_pack::graph_t>& clustered_index,
		const PairInfo<typename graph_pack::graph_t::EdgeId>& p_info) {
	int cnt = 0;
	auto pi_vector = clustered_index.GetEdgePairInfo(
			origin_gp.g.conjugate(p_info.second),
			origin_gp.g.conjugate(p_info.first));
	double tmpd = p_info.d + origin_gp.g.length(p_info.second)
			- origin_gp.g.length(p_info.first);
	for (auto pi_iter = pi_vector.begin(); pi_iter != pi_vector.end();
			++pi_iter) {
//					INFO("dist "<< tmpd << " versus "<<pi_iter->d<<" var "<<pi_iter->variance);
		if (abs(pi_iter->d - tmpd) < 0.01) {
			clustered_index.RemovePairInfo(*pi_iter);
			clustered_index.RemovePairInfo(BackwardInfo(*pi_iter));
			cnt += 2;
			INFO("Removed pi "<< origin_gp.g.int_id(pi_iter->first) << " "<< origin_gp.g.int_id(pi_iter->second)<<" d "<<pi_iter->d<<" var "<<pi_iter->variance);
		}
	}
	INFO("cnt += "<<cnt);
	return cnt;
}

template<class graph_pack>
int TreatPairPairInfo(const graph_pack& origin_gp,
		PairedInfoIndex<typename graph_pack::graph_t>& clustered_index,
		PairInfo<typename graph_pack::graph_t::EdgeId>& first_info,
		PairInfo<typename graph_pack::graph_t::EdgeId>& second_info,
		bool fill_missing, bool extensive_add = false) {

	size_t max_comparable_path = *cfg::get().ds.IS - K
			+ size_t(*cfg::get().ds.is_var);
	auto first_edge = first_info.second;
	auto first_weight = first_info.weight;
	if (fill_missing
			&& (first_info.d < -0.0001 || second_info.d < 0.0001
					|| first_info.d > second_info.d))
		return 0;
	auto second_edge = second_info.second;
	if (!fill_missing
			&& (first_info.d * second_info.d < 0.0001
					|| first_info.d > second_info.d))
		return 0;
	auto second_weight = second_info.weight;
	DEBUG(
			"Treating edges " << origin_gp.int_ids.ReturnIntId(first_edge)
					<< " " << origin_gp.int_ids.ReturnIntId(second_edge));
	vector < vector<typename graph_pack::graph_t::EdgeId> > paths;
	if (fill_missing) {
		paths = GetAllPathsBetweenEdges(
				origin_gp,
				first_edge,
				second_edge,
				llround(second_info.d - first_info.d)
						- origin_gp.g.length(first_edge) - 2,
				llround(second_info.d - first_info.d)
						- origin_gp.g.length(first_edge) + 2);
	} else {
		paths = GetAllPathsBetweenEdges(origin_gp, first_edge, second_edge, 0,
				*cfg::get().ds.IS - K + size_t(*cfg::get().ds.is_var));
	}
	vector < size_t > distances;
	for (auto paths_it = paths.begin(); paths_it != paths.end(); paths_it++) {
		distances.push_back(
				CummulativeLength<typename graph_pack::graph_t>(origin_gp.g,
						*paths_it));
	}
	//= callback.distances();
//	paths = GetAllPathsBetweenEdges(origin_gp, second_edge, first_edge);
//	for (auto paths_it = paths.begin(); paths_it != paths.end(); paths_it ++) {
//		distances.push_back(CummulativeLength<typename graph_pack::graph_t>(origin_gp.g, *paths_it));
//	}
	bool comparable = false;
	for (size_t i = 0; i < distances.size(); i++) {
		DEBUG(distances[i]);
		if (distances[i] < max_comparable_path) {
			comparable = true;
			break;
		}
	}
	if (!fill_missing) {
		if (comparable == false) {
			double ratio = (1.0 * second_weight) / first_weight;
			if (ratio > 1)
				ratio = 1 / ratio;
			if (first_weight > second_weight * 2) {
				clustered_index.RemovePairInfo(second_info);
				clustered_index.RemovePairInfo(BackwardInfo(second_info));
				DeleteConjugatePairInfo(origin_gp, clustered_index,
						second_info);
			} else if (second_weight > first_weight * 2) {
				clustered_index.RemovePairInfo(first_info);
				clustered_index.RemovePairInfo(BackwardInfo(first_info));
				DeleteConjugatePairInfo(origin_gp, clustered_index, first_info);
			}
			DEBUG(
					"contradictional paired info from edge "
							<< origin_gp.int_ids.ReturnIntId(first_info.first)
							<< " to edges "
							<< origin_gp.int_ids.ReturnIntId(first_edge)
							<< " and "
							<< origin_gp.int_ids.ReturnIntId(second_edge)
							<< "; weights ratio " << ratio);
			return 1;
		} else {
			DEBUG("no contradictions");
			return 0;
		}
	} else {
		if (paths.size() == 1 && first_info.variance == 0
				&& second_info.variance == 0
				&& (abs(
						distances[0] - second_info.d + first_info.d
								+ origin_gp.g.length(first_edge)) < 0.1)) {
			int added_info = 0;
			double tmpd = first_info.d + origin_gp.g.length(first_info.second);
			double w = (first_info.weight + second_info.weight) / 2;
			for (auto path_iter = paths.begin()->begin();
					path_iter != paths.begin()->end(); path_iter++) {
				if (TryToAddPairInfo(origin_gp, clustered_index,
						first_info.first, *path_iter, tmpd, w))
					added_info++;
				if (extensive_add) {
					if (TryToAddPairInfo(origin_gp, clustered_index,
							first_info.second, *path_iter, tmpd - first_info.d,
							w))
						added_info++;
					if (TryToAddPairInfo(origin_gp, clustered_index, *path_iter,
							second_info.second, second_info.d - tmpd, w))
						added_info++;
				}
				tmpd += origin_gp.g.length(*path_iter);
			}

			if (added_info) {

//				if (paths.begin()->size() != 0) {
				DEBUG("filled missing " << added_info << "edges");
				DEBUG(
						"while treating info from "
								<< origin_gp.int_ids.ReturnIntId(
										first_info.first) << " to "
								<< origin_gp.int_ids.ReturnIntId(first_edge)
								<< " "
								<< origin_gp.int_ids.ReturnIntId(second_edge));
//				}
				return added_info;
			}
		}
		return 0;
	}
}



template<class graph_pack>
bool isConsistent(const graph_pack& origin_gp, PairInfo<typename graph_pack::graph_t::EdgeId>& first_info,
		PairInfo<typename graph_pack::graph_t::EdgeId>& second_info) {
	if ((first_info.d < 0.0001 || second_info.d < 0.0001)
					|| (first_info.d > second_info.d)) return true;
//	size_t max_comparable_path = *cfg::get().ds.IS - K
//			+ size_t(*cfg::get().ds.is_var);
	int pi_distance = second_info.d - first_info.d;
	int first_length = origin_gp.g.length(first_info.second);
	double variance = first_info.variance + second_info.variance;

	auto first_edge = first_info.second;
	auto second_edge = second_info.second;
	INFO("   PI "<< origin_gp.g.int_id(first_info.first)<<" "<<origin_gp.g.int_id(first_info.second)<<" d "<<first_info.d<<"var "<<first_info.variance<<" tr "<< omp_get_thread_num());
	INFO("vs PI "<< origin_gp.g.int_id(second_info.first)<<" "<<origin_gp.g.int_id(second_info.second)<<" d "<<second_info.d<<"var "<<second_info.variance<<" tr "<< omp_get_thread_num());


	if (abs(pi_distance - first_length)<=variance){
		if (origin_gp.g.EdgeEnd(first_edge) == origin_gp.g.EdgeStart(second_edge)) return true;
		else {
			auto paths = GetAllPathsBetweenEdges(origin_gp, first_edge, second_edge, 0,
					pi_distance - first_length + variance);
			return (paths.size() > 0);
		}
	}
	else {
		if ((int)pi_distance > first_length){
			auto paths = GetAllPathsBetweenEdges(origin_gp, first_edge, second_edge, pi_distance - first_length - variance,
					pi_distance - first_length + variance);
			return (paths.size() > 0);
		}
		return false;
	}



//		vector < size_t > distances;
//	for (auto paths_it = paths.begin(); paths_it != paths.end(); paths_it++) {
//		distances.push_back(
//				CummulativeLength<typename graph_pack::graph_t>(origin_gp.g,
//						*paths_it));
//	}
//
//	TRACE("dist "<<distances<<"  tr"<< omp_get_thread_num());
//
//	bool comparable = false;
//	for (size_t i = 0; i < distances.size(); i++) {
//		TRACE(distances[i]);
//		if (distances[i] < pi_distance) {
//			comparable = true;
//			break;
//		}
//	}
//
//	return comparable;
//
}



template<class graph_pack>
void FindInconsistent(const graph_pack& origin_gp,
		PairedInfoIndex<typename graph_pack::graph_t>& clustered_index, std::vector<PairInfo<typename graph_pack::graph_t::EdgeId>>& edge_infos, PairInfoIndexData<typename graph_pack::graph_t::EdgeId>* pi) {

	for(size_t i = 0; i < edge_infos.size(); i++) {
		for(size_t j = 0; j < edge_infos.size(); j++) {
			if (i!=j){
				if (!isConsistent(origin_gp, edge_infos[i], edge_infos[j])){
					INFO("Inconsistent!!!")
					if (edge_infos[i].weight > edge_infos[j].weight) {
						pi->AddPairInfo(edge_infos[j]);
					} else  {
						pi->AddPairInfo(edge_infos[i]);
					}
				}
			}
		}
	}
}


template<class graph_pack>
int ParalelFillMissing(const graph_pack& origin_gp,
		PairedInfoIndex<typename graph_pack::graph_t>& clustered_index, size_t nthreads = 4){
	typedef typename graph_pack::graph_t Graph;
	int cnt = 0;
	INFO("Put infos to vector");
	std::vector<std::vector<PairInfo<typename graph_pack::graph_t::EdgeId>>> infos;
	for (auto e_iter = origin_gp.g.SmartEdgeBegin(); !e_iter.IsEnd();
			++e_iter) {
//		if (origin_gp.g.length(*e_iter)	>= cfg::get().rr.max_repeat_length)
			infos.push_back(clustered_index.GetEdgeInfo(*e_iter));
	}

	std::vector< PairedInfoIndex<Graph>* > to_add(nthreads);
    for (size_t i = 0; i < nthreads; ++i) {
        to_add[i] = new PairedInfoIndex<Graph>(origin_gp.g);
    }
    SplitPathConstructor<Graph> spc(origin_gp.g);
    INFO("Start threads");
	size_t n = 0;
	#pragma omp parallel num_threads(nthreads)
	{
		#pragma omp for
		for (size_t i = 0; i < infos.size(); ++i)
		{
//			FindInconsistent<graph_pack>(origin_gp, clustered_index, infos[i], to_add[omp_get_thread_num()]);
			vector<PathInfoClass<Graph>> paths = spc.ConvertEdgePairInfoToSplitPathes(infos[i]);
			for (auto iter = paths.begin(); iter != paths.end(); iter ++) {
				INFO("Path "<<iter->PrintPath(origin_gp.g));
				for (auto pi_iter = iter->begin(); pi_iter != iter->end(); pi_iter ++) {
					TryToAddPairInfo(origin_gp, *to_add[omp_get_thread_num()], *pi_iter);
				}
			}
			VERBOSE_POWER_T(++n, 100, " edge pairs processed. Cur thread "<<omp_get_thread_num());

		}
	}

	INFO("Threads finished");

	for (size_t i = 0; i < nthreads; ++i) {
		for (auto add_iter = (*to_add[i]).begin(); add_iter != (*to_add[i]).end(); ++add_iter) {
			auto pinfos = *add_iter;
			for (auto add_pi_iter = pinfos.begin(); add_pi_iter != pinfos.end(); ++add_pi_iter) {
				if (TryToAddPairInfo(origin_gp, clustered_index, *add_pi_iter))
				{
					INFO("PI added "<<cnt);
					cnt++;
				}
			}
		}
        delete to_add[i];
	}
	INFO("Clean finished");
	return cnt;
}


template<class graph_pack>
int ParalelRemoveContraditional(const graph_pack& origin_gp,
		PairedInfoIndex<typename graph_pack::graph_t>& clustered_index, size_t nthreads = 4){

	int cnt = 0;
	INFO("Put infos to vector");
	std::vector<std::vector<PairInfo<typename graph_pack::graph_t::EdgeId>>> infos;
	for (auto e_iter = origin_gp.g.SmartEdgeBegin(); !e_iter.IsEnd();
			++e_iter) {
		if (origin_gp.g.length(*e_iter)	>= cfg::get().rr.max_repeat_length)
			infos.push_back(clustered_index.GetEdgeInfo(*e_iter));
	}

	std::vector< PairInfoIndexData<typename graph_pack::graph_t::EdgeId>* > to_remove(nthreads);
    for (size_t i = 0; i < nthreads; ++i) {
        to_remove[i] = new PairInfoIndexData<typename graph_pack::graph_t::EdgeId>();
    }
	INFO("Start threads");
	size_t n = 0;
	#pragma omp parallel num_threads(nthreads)
	{
		#pragma omp for
		for (size_t i = 0; i < infos.size(); ++i)
		{
			FindInconsistent<graph_pack>(origin_gp, clustered_index, infos[i], to_remove[omp_get_thread_num()]);
            VERBOSE_POWER_T(++n, 100, " edge pairs processed. Cur thread "<<omp_get_thread_num());

		}
	}

	INFO("Threads finished");

	for (size_t i = 0; i < nthreads; ++i) {
		for (auto remove_iter = (*to_remove[i]).begin(); remove_iter != (*to_remove[i]).end(); ++remove_iter) {
			cnt += DeleteIfExist<graph_pack>(origin_gp, clustered_index, *remove_iter);
			cnt += DeleteConjugatePairInfo<graph_pack>(origin_gp, clustered_index, *remove_iter);
		}
        delete to_remove[i];
	}
	INFO("Clean finished");
	return cnt;
}


template<class graph_pack>
void ParalelCorrectPairedInfo(const graph_pack& origin_gp,
		PairedInfoIndex<typename graph_pack::graph_t>& clustered_index, size_t nthreads = 4) {
	int missing_paired_info_count = 0;
	int extra_paired_info_count = 0;

	extra_paired_info_count = ParalelRemoveContraditional(origin_gp, clustered_index, nthreads);
	missing_paired_info_count = ParalelFillMissing(origin_gp, clustered_index, nthreads);

	INFO(
			"Paired info stats: missing = " << missing_paired_info_count
					<< "; contradictional = " << extra_paired_info_count);
}



template<class graph_pack>
void CorrectPairedInfo(const graph_pack& origin_gp,
		PairedInfoIndex<typename graph_pack::graph_t>& clustered_index,
		bool clean = true, bool add = true) {
	size_t k = graph_pack::k_value;
	size_t delta = size_t(*cfg::get().ds.is_var);
	size_t max_comparable_path = *cfg::get().ds.IS + delta - k;
	int missing_paired_info_count = 0;
	int extra_paired_info_count = 0;
	int long_edges_count = 0;

	DEBUG("Using max path cutoff = " << max_comparable_path);
	for (auto e_iter = origin_gp.g.SmartEdgeBegin(); !e_iter.IsEnd();
			++e_iter) {
		if (origin_gp.g.length(*e_iter) >= cfg::get().rr.max_repeat_length)
			long_edges_count++;
		auto pi = clustered_index.GetEdgeInfo(*e_iter);
		for (auto i_iter = pi.begin(); i_iter != pi.end(); ++i_iter) {
			for (auto j_iter = pi.begin(); j_iter != pi.end(); ++j_iter) {
				if (i_iter != j_iter) {
					PairInfo<typename graph_pack::graph_t::EdgeId> first_info =
							*i_iter;
					PairInfo<typename graph_pack::graph_t::EdgeId> second_info =
							*j_iter;
					if (origin_gp.g.length(*e_iter)
							>= /* *cfg::get().ds.RL * */2 && add) { //TODO: change to something reasonable.
						missing_paired_info_count += TreatPairPairInfo<
								graph_pack>(origin_gp, clustered_index,
								first_info, second_info, 1);
					}
					if (origin_gp.g.length(*e_iter)
							>= cfg::get().rr.max_repeat_length && clean) {
						extra_paired_info_count +=
								TreatPairPairInfo<graph_pack>(origin_gp,
										clustered_index, first_info,
										second_info, 0);
					}
				}
			}
		}
	}
	INFO(
			"Paired info stats: missing = " << missing_paired_info_count
					<< "; contradictional = " << extra_paired_info_count);
}
template<class graph_pack>
void process_resolve_repeats(graph_pack& origin_gp,
		PairedInfoIndex<typename graph_pack::graph_t>& clustered_index,
		graph_pack& resolved_gp, const string& graph_name,
		const string& subfolder = "", bool output_contigs = true) {

	EdgeLabelHandler<typename graph_pack::graph_t> labels_after(resolved_gp.g,
			origin_gp.g);

	process_resolve_repeats(origin_gp, clustered_index, resolved_gp, graph_name,
			labels_after, subfolder, output_contigs);
}

template<class graph_pack>
void component_statistics(graph_pack & conj_gp, int component_id,
		PairedInfoIndex<typename graph_pack::graph_t>& clustered_index) {

	string graph_name = ConstructComponentName("graph_", component_id).c_str();
	string component_name = cfg::get().output_dir + "graph_components/"
			+ graph_name;
	//component output
	string table_name = cfg::get().output_dir + "graph_components/tables/";
	make_dir(table_name);
	table_name += graph_name;
	set<typename graph_pack::graph_t::EdgeId> incoming_edges;
	set<typename graph_pack::graph_t::EdgeId> outgoing_edges;
	for (auto iter = conj_gp.g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
		typename graph_pack::graph_t::VertexId start = conj_gp.g.EdgeStart(
				*iter);
		typename graph_pack::graph_t::VertexId end = conj_gp.g.EdgeEnd(*iter);
		if (conj_gp.g.length(*iter) > *cfg::get().ds.IS + 100) {

			if (conj_gp.g.IsDeadStart(
					start) /*&& conj_gp.g.CheckUniqueOutgoingEdge(start)*/) {
				incoming_edges.insert(*iter);
			} else if (conj_gp.g.IsDeadEnd(
					end)/* && conj_gp.g.CheckUniqueIncomingEdge(end)*/) {
				outgoing_edges.insert(*iter);
			} else {
				WARN(
						"strange long edge in component " << component_name
								<< " , edge_id "
								<< conj_gp.int_ids.ReturnIntId(*iter));
			}
		}
	}
	INFO("incoming- outgoint set formed");
	int flag = 1;
	for (auto inc_iter = incoming_edges.begin();
			inc_iter != incoming_edges.end(); ++inc_iter) {
		int count = 0;
		for (auto out_iter = outgoing_edges.begin();
				out_iter != outgoing_edges.end(); ++out_iter) {
			if (clustered_index.GetEdgePairInfo(*inc_iter, *out_iter).size()
					== 1)
				count++;
		}
		if (count != 1)
			flag = 0;
	}
	FILE* file;
	if (flag)
		file = fopen((table_name + ".tbl_good").c_str(), "w");
	else
		file = fopen((table_name + ".tbl").c_str(), "w");

	INFO("Saving in-out table , " << component_name << " created");
	VERIFY(file != NULL);
	fprintf(file, "%7c", ' ');

	for (auto out_iter = outgoing_edges.begin();
			out_iter != outgoing_edges.end(); ++out_iter)
		fprintf(file, " %7d", conj_gp.int_ids.ReturnIntId(*out_iter));
	fprintf(file, "\n");

	for (auto inc_iter = incoming_edges.begin();
			inc_iter != incoming_edges.end(); ++inc_iter) {
		fprintf(file, " %7d", conj_gp.int_ids.ReturnIntId(*inc_iter));
		for (auto out_iter = outgoing_edges.begin();
				out_iter != outgoing_edges.end(); ++out_iter) {
			char c;
			if (clustered_index.GetEdgePairInfo(*inc_iter, *out_iter).size()
					== 0)
				c = '0';
			else
				c = 'X';
			fprintf(file, "%7c", c);
		}
		fprintf(file, "\n");
	}

	fprintf(file, "\n");
	for (auto inc_iter = incoming_edges.begin();
			inc_iter != incoming_edges.end(); ++inc_iter)
		fprintf(file, " %7d", conj_gp.int_ids.ReturnIntId(*inc_iter));
	fprintf(file, "\n");
	for (auto out_iter = outgoing_edges.begin();
			out_iter != outgoing_edges.end(); ++out_iter)
		fprintf(file, " %7d", conj_gp.int_ids.ReturnIntId(*out_iter));
	fprintf(file, "\n");

	fclose(file);

}

void resolve_conjugate_component(int component_id, const Sequence& genome) {
	conj_graph_pack conj_gp(genome, cfg::get().pos.max_single_gap,
			cfg::get().pos.careful_labeling);
	paired_info_index paired_index(conj_gp.g/*, 5.*/);
	paired_info_index clustered_index(conj_gp.g);

	INFO("Resolve component " << component_id);

	string graph_name = ConstructComponentName("graph_", component_id).c_str();
	string component_name = cfg::get().output_dir + "graph_components/"
			+ graph_name;

	ScanWithClusteredIndex(component_name, conj_gp, clustered_index);

	component_statistics(conj_gp, component_id, clustered_index);

	conj_graph_pack resolved_gp(genome, cfg::get().pos.max_single_gap,
			cfg::get().pos.careful_labeling);
	string sub_dir = "resolve_components/";

	string resolved_name = cfg::get().output_dir + "resolve_components"
			+ "/resolve_" + graph_name + "/";
	make_dir(resolved_name);
	process_resolve_repeats(conj_gp, clustered_index, resolved_gp, graph_name,
			sub_dir, false);
}

void resolve_nonconjugate_component(int component_id, const Sequence& genome) {
//	nonconj_graph_pack nonconj_gp(genome);
//	PairedInfoIndex<nonconj_graph_pack::graph_t> clustered_index(nonconj_gp.g);
//
//	INFO("Resolve component "<<component_id);
//
//	string graph_name = ConstructComponentName("graph_", component_id).c_str();
//	string component_name = cfg::get().output_dir + "graph_components/"
//			+ graph_name;
//
//	ScanWithClusteredIndex(component_name, nonconj_gp, clustered_index);
//
//	component_statistics(nonconj_gp, component_id, clustered_index);
//
//	nonconj_graph_pack resolved_gp(genome);
//	string sub_dir = "resolve_components/";
//
//	string resolved_name = cfg::get().output_dir + "resolve_components"
//			+ "/resolve_" + graph_name + "/";
//	make_dir(resolved_name);
//	process_resolve_repeats(nonconj_gp, clustered_index, resolved_gp,
//			graph_name, sub_dir, false);
}

void resolve_with_jumps(conj_graph_pack& gp, PairedInfoIndex<Graph>& index,
		const paired_info_index& jump_index) {
	VERIFY(cfg::get().andrey_params.write_contigs);
	resolve_repeats_ml(gp, index, cfg::get().output_dir + "jump_resolve/",
			cfg::get().andrey_params,
			boost::optional<const paired_info_index&>(jump_index));
}

void prepare_jump_index(const Graph& g, const paired_info_index& raw_jump_index,
		paired_info_index& jump_index) {
	JumpingEstimator<Graph> estimator(raw_jump_index);
	paired_info_index clustered_jump_index(g);
	estimator.Estimate(clustered_jump_index);

	JumpingNormilizerFunction<Graph> nf(g, *cfg::get().ds.RL, 500);
	PairedInfoNormalizer<Graph> normalizer(clustered_jump_index, nf);
	paired_info_index normalized_jump_index(g);
	normalizer.FillNormalizedIndex(normalized_jump_index);

	JumpingPairInfoChecker<Graph> filter(g, 300, 100, 100);
	filter.Filter(normalized_jump_index, jump_index);
}

void resolve_repeats() {
	Sequence genome = cfg::get().ds.reference_genome;

	conj_graph_pack conj_gp(genome, cfg::get().pos.max_single_gap,
			cfg::get().pos.careful_labeling);
	paired_info_index paired_index(conj_gp.g, cfg::get().online_clust_rad);
	paired_info_index clustered_index(conj_gp.g);

	exec_distance_estimation(conj_gp, paired_index, clustered_index);

	DEBUG("Online clusterization rad = " << cfg::get().online_clust_rad);
	if (cfg::get().pos.late_threading) {
		FillPos(conj_gp, conj_gp.genome, "10");
		FillPos(conj_gp, !conj_gp.genome, "11");
		if (!cfg::get().pos.contigs_for_threading.empty()
				&& fileExists(cfg::get().pos.contigs_for_threading)) {
			FillPos(conj_gp, cfg::get().pos.contigs_for_threading, 10000);
		}
	}

	//todo refactor labeler creation
	total_labeler_graph_struct graph_struct(conj_gp.g, &conj_gp.int_ids,
			&conj_gp.edge_pos);
	total_labeler tot_lab(&graph_struct);
	EdgeQuality<Graph> quality_labeler(conj_gp.g, conj_gp.index,
			conj_gp.kmer_mapper, conj_gp.genome);
//	OutputWrongContigs<K>(conj_gp, 1000, "contamination.fasta");
	CompositeLabeler<Graph> labeler(tot_lab, quality_labeler);
	detail_info_printer printer(conj_gp, labeler, cfg::get().output_dir,
			"graph.dot");
	printer(ipp_before_repeat_resolution);

	if (!cfg::get().paired_mode
			|| cfg::get().rm == debruijn_graph::resolving_mode::rm_none) {
		OutputContigs(conj_gp.g, cfg::get().output_dir + "final_contigs.fasta");
		return;
	}

	//tSeparatedStats(conj_gp, conj_gp.genome, clustered_index);

	SAM_before_resolve(conj_gp);

	INFO("STAGE == Resolving Repeats");

	if (cfg::get().rm == debruijn_graph::resolving_mode::rm_split) {
		int number_of_components = 0;

		if (cfg::get().componential_resolve) {
			make_dir(cfg::get().output_dir + "graph_components" + "/");
			number_of_components = PrintGraphComponents(
					cfg::get().output_dir + "graph_components/graph_", conj_gp,
					*cfg::get().ds.IS + 100, clustered_index);
			INFO("number of components " << number_of_components);
		}

		if (cfg::get().rr.symmetric_resolve) {
			if (cfg::get().use_multithreading) {
				ParalelCorrectPairedInfo(conj_gp, clustered_index, cfg::get().max_threads);
				ParalelCorrectPairedInfo(conj_gp, clustered_index, cfg::get().max_threads);
			} else {
				CorrectPairedInfo(conj_gp, clustered_index, true, false);
//				CorrectPairedInfo(conj_gp, clustered_index, true, false);
//				CorrectPairedInfo(conj_gp, clustered_index, true, false);
//				CorrectPairedInfo(conj_gp, clustered_index);
//				CorrectPairedInfo(conj_gp, clustered_index);
//				CorrectPairedInfo(conj_gp, clustered_index);
//				CorrectPairedInfo(conj_gp, clustered_index);
			}

			save_distance_filling(conj_gp, clustered_index, clustered_index);

			conj_graph_pack resolved_gp(genome, cfg::get().pos.max_single_gap,
					cfg::get().pos.careful_labeling);

			EdgeLabelHandler<conj_graph_pack::graph_t> labels_after(
					resolved_gp.g, conj_gp.g);

			process_resolve_repeats(conj_gp, clustered_index, resolved_gp,
					"graph", labels_after, "", true);

			SAM_after_resolve(conj_gp, resolved_gp, labels_after);

			if (cfg::get().componential_resolve) {
				make_dir(cfg::get().output_dir + "resolve_components" + "/");
				for (int i = 0; i < number_of_components; i++) {
					resolve_conjugate_component(i + 1, genome);
				}
			}
		} else {
//			nonconj_graph_pack origin_gp(conj_gp.genome);
//			PairedInfoIndex<nonconj_graph_pack::graph_t> orig_clustered_idx(
//					origin_gp.g);
//			Convert(conj_gp, clustered_index, origin_gp, orig_clustered_idx);
//			nonconj_graph_pack resolved_gp(conj_gp.genome);
//			process_resolve_repeats(origin_gp, orig_clustered_idx, resolved_gp,
//					"graph");
//			if (cfg::get().componential_resolve) {
//				make_dir(cfg::get().output_dir + "resolve_components" + "/");
//				for (int i = 0; i < number_of_components; i++) {
//					resolve_nonconjugate_component(i + 1, genome);
//				}
//			}
		}
	}

	//todo magic constants!!!
	if (cfg::get().rm == debruijn_graph::resolving_mode::rm_jump) {
		if (!cfg::get().jump.load) {
			INFO("Going to count jumping paired info");
			VERIFY(
					cfg::get().ds.jumping_first && cfg::get().ds.jumping_second
							&& cfg::get().ds.jump_is);
			checkFileExistenceFATAL(
					cfg::get().input_dir + (*cfg::get().ds.jumping_first));
			checkFileExistenceFATAL(
					cfg::get().input_dir + (*cfg::get().ds.jumping_second));
			paired_info_index raw_jump_index(conj_gp.g, 1000.);
			io::PairedEasyReader jump_stream(
					make_pair(
							cfg::get().input_dir
									+ (*cfg::get().ds.jumping_first),
							cfg::get().input_dir
									+ (*cfg::get().ds.jumping_second)),
					*cfg::get().ds.jump_is, true);

//			cout << "eof " << jump_stream.eof() << endl;
//			cout << "START HERE" << endl;
//			io::PairedRead paired_read;
//			while (!jump_stream.eof()) {
//				jump_stream >> paired_read;
//				cout << "HERE " << paired_read.first().sequence() << endl;
//			}
//			cout << "END HERE" << endl;

			io::ISCorruptingWrapper wrapped_jump_stream(jump_stream, 1e6);
			std::vector<io::IReader<io::PairedRead>*> stream;
			stream.push_back(&wrapped_jump_stream);

			FillPairedIndexWithReadCountMetric<K>(conj_gp.g, conj_gp.int_ids,
					conj_gp.index, conj_gp.kmer_mapper, raw_jump_index, stream);
//			FillPairedIndex<K>(conj_gp.g,
//					conj_gp.index, raw_jump_index,
//					wrapped_jump_stream);

			ConjugateDataPrinter<Graph> printer(conj_gp.g, conj_gp.int_ids);
			printer.savePaired(cfg::get().output_dir + "jump_raw",
					raw_jump_index);

			paired_info_index jump_index(conj_gp.g);
			prepare_jump_index(conj_gp.g, raw_jump_index, jump_index);

			printer.savePaired(cfg::get().output_dir + "jump_cleared",
					jump_index);
			resolve_with_jumps(conj_gp, clustered_index, jump_index);
		} else {
			ConjugateDataScanner<Graph> scanner(conj_gp.g, conj_gp.int_ids);
			paired_info_index jump_index(conj_gp.g);
			scanner.loadPaired(cfg::get().output_dir + "../jump_cleared",
					jump_index);
			resolve_with_jumps(conj_gp, clustered_index, jump_index);

//			ConjugateDataScanner<Graph> scanner(conj_gp.g, conj_gp.int_ids);
//			paired_info_index raw_jump_index(conj_gp.g);
//			scanner.loadPaired(cfg::get().output_dir + "../jump_raw",
//					raw_jump_index);
//
//			paired_info_index jump_index(conj_gp.g);
//			prepare_jump_index(conj_gp.g, raw_jump_index, jump_index);
//
//			ConjugateDataPrinter<Graph> printer(conj_gp.g, conj_gp.int_ids);
//			printer.savePaired(cfg::get().output_dir + "jump_cleared",
//					jump_index);
//
//			resolve_with_jumps(conj_gp, clustered_index, jump_index);
		}
	}

	if (cfg::get().rm == debruijn_graph::resolving_mode::rm_path_extend) {
		resolve_repeats_ml(conj_gp, clustered_index,
				cfg::get().output_dir + "alt_resolve/",
				cfg::get().andrey_params);
	}

	if (cfg::get().rm == debruijn_graph::resolving_mode::rm_combined
			&& !cfg::get().path_set_graph) {
		INFO("Combined resolving started");

		std::string graph_name = "resolved_graph";

		conj_graph_pack resolved_gp(genome, cfg::get().pos.max_single_gap,
				cfg::get().pos.careful_labeling);

		EdgeLabelHandler<conj_graph_pack::graph_t> labels_after(resolved_gp.g,
				conj_gp.g);

		process_resolve_repeats(conj_gp, clustered_index, resolved_gp, "graph",
				labels_after);

		PairedInfoIndex<conj_graph_pack::graph_t> resolved_graph_paired_info(
				resolved_gp.g);
		ProduceResolvedPairedInfo(conj_gp, clustered_index, resolved_gp,
				labels_after, resolved_graph_paired_info);
		SaveResolvedPairedInfo(resolved_gp, resolved_graph_paired_info,
				"resolved", "saves/");

		resolve_repeats_ml(resolved_gp, resolved_graph_paired_info,
				cfg::get().output_dir + "combined_resolve/",
				cfg::get().andrey_params);

		INFO("Combined resolving finished");
	}

}

void exec_repeat_resolving() {
	if (cfg::get().entry_point <= ws_repeats_resolving) {
		resolve_repeats();
		//todo why nothing to save???
		// nothing to save yet
	} else {
		INFO("Loading Repeat Resolving");
		INFO("Nothing to load");
		// nothing to load
	}
}

} // debruijn_graph

