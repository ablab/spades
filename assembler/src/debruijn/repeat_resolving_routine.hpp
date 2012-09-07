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

#include "internal_aligner.hpp"
#include "omni/loop_killer.hpp"
#include "path_utils.hpp"
#include "pair_info_improver.hpp"

#include "path_extend/path_extend_launch.hpp"

typedef io::CarefulFilteringReaderWrapper<io::SingleRead> CarefulFilteringStream;

namespace debruijn_graph {

void resolve_repeats(PairedReadStream& stream, const Sequence& genome);
} // debruijn_graph

// TODO move impl to *.cpp

namespace debruijn_graph {

template<class gp_t>
void WriteGraphPack(gp_t& gp, const string& file_name) {
	ofstream filestr(file_name);
	CompositeGraphColorer<typename gp_t::graph_t> colorer(
			new FixedColorer<typename gp_t::graph_t::VertexId>("white"),
			new PositionsEdgeColorer<typename gp_t::graph_t>(gp.g, gp.edge_pos));

	EdgeQuality<typename gp_t::graph_t> edge_qual(gp.g, gp.index, gp.kmer_mapper, gp.genome);
	total_labeler_graph_struct graph_struct(gp.g, &gp.int_ids, &gp.edge_pos);
	total_labeler tot_lab(&graph_struct);
	CompositeLabeler<Graph> labeler(tot_lab, edge_qual);

	DotGraphPrinter<typename gp_t::graph_t> g_print(gp.g, labeler, colorer, " ", filestr);
	SimpleGraphVisualizer<typename gp_t::graph_t> gv(gp.g, g_print);
	gv.Visualize();
}

void save_distance_filling(conj_graph_pack& gp, paired_info_index& paired_index,
		paired_info_index& clustered_index) {
	if (cfg::get().make_saves) {
        string p = path::append_path(cfg::get().output_saves, "distance_filling");
        PrintAll(p, gp, paired_index, clustered_index);
        write_estimated_params(p);
	}
}

bool try_load_distance_filling(conj_graph_pack& gp, paired_info_index& clustered_index,
        path::files_t* used_files) {

    string p = path::append_path(cfg::get().load_from, "distance_filling");

    FILE* file = fopen((p + ".grp").c_str(), "r");
    if (file == NULL) {
        return false;
    }
    fclose(file);

    used_files->push_back(p);

    clustered_index.Clear();
    ScannerTraits<conj_graph_pack::graph_t>::Scanner scanner(gp.g,
                gp.int_ids);
    ScanClusteredIndex(p, scanner, clustered_index);

    return true;
}


void distance_filling(conj_graph_pack& gp, paired_info_index& paired_index,
        paired_info_index& clustered_index) {

    path::files_t used_files;
    if (try_load_distance_filling(gp, clustered_index, &used_files)) {

        link_files_by_prefix(used_files, cfg::get().output_saves);
        INFO("Distance filling saves detected and loaded");
    }
    else {
        INFO("Filling paired information");

        PairInfoInprover<conj_graph_pack::graph_t> pi_imp(gp.g);
        pi_imp.ImprovePairedInfo(clustered_index,
                    cfg::get().use_multithreading, cfg::get().max_threads);

        save_distance_filling(gp, paired_index, clustered_index);
    }
}


void save_resolved(conj_graph_pack& resolved_gp,
        paired_info_index& resolved_graph_paired_info,
        paired_info_index& resolved_graph_paired_info_cl) {

    if (cfg::get().make_saves) {
        string p = path::append_path(cfg::get().output_saves, "split_resovled");
        PrintAll(p, resolved_gp, resolved_graph_paired_info, resolved_graph_paired_info_cl);
        write_estimated_params(p);
    }
}

template<class graph_pack>
void SelectReadsForConsensus(graph_pack& etalon_gp,
		typename graph_pack::graph_t& cur_graph,
		EdgeLabelHandler<typename graph_pack::graph_t>& LabelsAfter,
		const EdgeIndex<typename graph_pack::graph_t>& index,
		vector<ReadStream *>& reads, string& consensus_output_dir, size_t k) {
	INFO("ReadMapping started");
	map<typename graph_pack::graph_t::EdgeId, int> contigNumbers;
	int cur_num = 0;
	FillContigNumbers(contigNumbers, cur_graph);
	for (auto iter = etalon_gp.g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
		DEBUG(
				"Edge number:" << etalon_gp.int_ids.ReturnIntId(*iter) << " is contained in contigs");
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
		SingleReadMapper<typename graph_pack::graph_t> rm(etalon_gp.g, index,
				k);
		DEBUG("mapping reads from pair " << i);
		while (!reads[i - 1]->eof()) {
			io::SingleRead cur_read;

			(*reads[i - 1]) >> cur_read;
			vector<typename graph_pack::graph_t::EdgeId> res =
					rm.GetContainingEdges(cur_read);
			read_num++;
			TRACE(
					read_num << " mapped to" << res.size() << " contigs :, read" << cur_read.sequence());
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
								"No edges containing" << etalon_gp.int_ids.ReturnIntId( res[ii]));
			}
		}
	}
}


void SAMAfterResolve(conj_graph_pack& conj_gp, conj_graph_pack& resolved_gp,
		EdgeLabelHandler<conj_graph_pack::graph_t> &labels_after) {


	io::OffsetType offset_type = EvaluateOffset();
	string OutputFileName = (cfg::get().run_mode) ? cfg::get().output_dir + "align_after_RR.sam": cfg::get().output_base + "contigs.sam";


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
			typedef NewExtendedSequenceMapper<Graph> SequenceMapper;
			SequenceMapper mapper(conj_gp.g, conj_gp.index, conj_gp.kmer_mapper,
					conj_gp.k_value + 1);

			bool print_quality = (
					cfg::get().sw.print_quality ?
							*cfg::get().sw.print_quality : false);
			OriginalReadsResolvedInternalAligner<ConjugateDeBruijnGraph,
					SequenceMapper> Aligner(resolved_gp.k_value, resolved_gp.g,
					conj_gp.g, mapper, labels_after, cfg::get().sw.adjust_align,
					cfg::get().sw.output_map_format,
					cfg::get().sw.output_broken_pairs, print_quality);
			Aligner.AlignPairedReads(*original_paired_reads, *paired_reads,
					OutputFileName);
		}
	} else {
		auto paired_reads = paired_easy_reader(false, 0, false, false, false,
				offset_type);
		auto single_reads = single_easy_reader(false, false, false,
				offset_type);

		typedef NewExtendedSequenceMapper<Graph> SequenceMapper;
		SequenceMapper mapper(conj_gp.g, conj_gp.index, conj_gp.kmer_mapper,
				conj_gp.k_value + 1);

		bool print_quality = (
				cfg::get().sw.print_quality ?
						*cfg::get().sw.print_quality : false);
		ResolvedInternalAligner<ConjugateDeBruijnGraph, SequenceMapper> Aligner(
				resolved_gp.k_value, resolved_gp.g, conj_gp.g, mapper,
				labels_after, cfg::get().sw.adjust_align,
				cfg::get().sw.output_map_format,
				cfg::get().sw.output_broken_pairs, print_quality);
		if (cfg::get().sw.align_only_paired)
			Aligner.AlignPairedReads(*paired_reads,
					OutputFileName);
		else
			Aligner.AlignReads(*paired_reads, *single_reads,
					OutputFileName);

	}
}

template<class graph_pack>
void CleanIsolated(graph_pack& gp) {
	for (auto iter = gp.g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
		typename graph_pack::graph_t::VertexId start, end;
		start = gp.g.EdgeStart(*iter);
		end = gp.g.EdgeEnd(*iter);
		TRACE(
				gp.g.CheckUniqueOutgoingEdge(start) << " " << gp.g.IsDeadStart(start) << " " << gp.g.CheckUniqueIncomingEdge(end) << " " << gp.g.IsDeadEnd(end));
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

	s += debruijn_config::estimation_mode_name(cfg::get().est_mode) + "_est_";

	s += "k";
	s += ToString(cfg::get().K);
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
void ProduceResolvedPairedInfo(graph_pack& origin_gp,
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
void SaveResolvedPairedInfo(graph_pack& resolved_gp,
		PairedInfoIndex<typename graph_pack::graph_t> resolved_graph_paired_info,
		const string& graph_name, const string& subfolder) {
	if (cfg::get().make_saves) {
		std::string rr_filename;
		if (subfolder.size()) {
			INFO("Saving graph and paired info to subfolder " << subfolder);
			rr_filename = (cfg::get().output_dir + subfolder) + graph_name;
		} else {
			INFO("Saving graph and paired info");
            string p = path::append_path(cfg::get().output_saves, graph_name);
            rr_filename = p;
		}
        PrintWithClusteredIndex(rr_filename, resolved_gp, resolved_graph_paired_info);
		DEBUG("Saved");
	}
}

template<class graph_pack>
void process_resolve_repeats(graph_pack& origin_gp,
		PairedInfoIndex<typename graph_pack::graph_t>& clustered_index,
		graph_pack& resolved_gp, const string& graph_name,
		EdgeLabelHandler<typename graph_pack::graph_t>& labels_after,
		const string& subfolder = "", bool output_contigs = true, bool kill_loops = true) {

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
		
	//EdgeQuality<Graph> quality_lab_(origin_gp.g, origin_gp.index,
			//origin_gp.kmer_mapper, origin_gp.genome);
	//CompositeLabeler<Graph> lab_(tot_labeler_before, quality_lab_);

    //omnigraph::WriteSimple(
				//origin_gp.g,
				//lab_,
				//cfg::get().output_dir + subfolder + graph_name
						//+ "_2_before.dot", "no_repeat_graph");

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
		ResolveRepeats(origin_gp.g, origin_gp.int_ids, clustered_index,
				origin_gp.edge_pos, resolved_gp.g, resolved_gp.int_ids,
				resolved_gp.edge_pos,
				cfg::get().output_dir + subfolder + "resolve_" + graph_name
						+ "/", labels_after, cfg::get().developer_mode);

		//Generating paired info for resolved graph
//		PairedInfoIndex<typename graph_pack::graph_t> resolved_graph_paired_info(resolved_gp.g);
//		ProduceResolvedPairedInfo(origin_gp, clustered_index, resolved_gp, labels_after, resolved_graph_paired_info);
//		SaveResolvedPairedInfo(resolved_gp, resolved_graph_paired_info, graph_name + "_resolved", subfolder);
		//Paired info for resolved graph generated

	}

	if (cfg::get().output_nonfinal_contigs && output_contigs) {
		OutputContigs(resolved_gp.g,
				cfg::get().output_dir + "after_rr_before_simplify" + postfix);
	}
	INFO("Running total labeler");

	total_labeler_gs graph_struct_after(resolved_gp.g, &resolved_gp.int_ids,
			&resolved_gp.edge_pos, &labels_after);
	total_labeler tot_labeler_after(&graph_struct_after, &graph_struct_before);
	if (cfg::get().output_pictures) {
		omnigraph::WriteSimple(resolved_gp.g, tot_labeler_after,
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
				"Tip clipping iteration " << i << " (0-indexed) out of " << iters << ":");

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

//	OnlineVisualizer online(resolved_gp);
//	online.run();
	if (kill_loops) {
		SimpleLoopKiller<typename graph_pack::graph_t> lk(resolved_gp.g,
				cfg::get().rr.max_repeat_length, 6);
		lk.KillAllLoops();
	}

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
		OutputContigs(resolved_gp.g,
				cfg::get().output_dir + "final_contigs.fasta");
		cfg::get_writable().final_contigs_file = cfg::get().output_dir
				+ "final_contigs.fasta";
	}

	if (cfg::get().output_pictures) {
		omnigraph::WriteSimple(resolved_gp.g, tot_labeler_after,
				cfg::get().output_dir + subfolder + graph_name
						+ "_4_cleared.dot", "no_repeat_graph");
		ofstream filestr(
				cfg::get().output_dir + subfolder + graph_name
						+ "_4_cleared_colored.dot");
		CompositeGraphColorer<typename graph_pack::graph_t> colorer(
				new FixedColorer<typename graph_pack::graph_t::VertexId>(
						"white"),
				new PositionsEdgeColorer<typename graph_pack::graph_t>(
						resolved_gp.g, resolved_gp.edge_pos));
		DotGraphPrinter<typename graph_pack::graph_t> gp(resolved_gp.g,
				tot_labeler_after, colorer, " ", filestr);
		SimpleGraphVisualizer<typename graph_pack::graph_t> gv(resolved_gp.g,
				gp);
		gv.Visualize();
		filestr.close();

	}
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
						"strange long edge in component " << component_name << " , edge_id " << conj_gp.int_ids.ReturnIntId(*iter));
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
		fprintf(file, " %7zu", conj_gp.int_ids.ReturnIntId(*out_iter));
	fprintf(file, "\n");

	for (auto inc_iter = incoming_edges.begin();
			inc_iter != incoming_edges.end(); ++inc_iter) {
		fprintf(file, " %7zu", conj_gp.int_ids.ReturnIntId(*inc_iter));
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
		fprintf(file, " %7zu", conj_gp.int_ids.ReturnIntId(*inc_iter));
	fprintf(file, "\n");
	for (auto out_iter = outgoing_edges.begin();
			out_iter != outgoing_edges.end(); ++out_iter)
		fprintf(file, " %7zu", conj_gp.int_ids.ReturnIntId(*out_iter));
	fprintf(file, "\n");

	fclose(file);

}

void resolve_conjugate_component(int component_id, const Sequence& genome) {
	conj_graph_pack conj_gp(cfg::get().K,
                          cfg::get().output_dir,
                          genome, cfg::get().pos.max_single_gap,
                          cfg::get().pos.careful_labeling);
	paired_info_index paired_index(conj_gp.g/*, 5.*/);
	paired_info_index clustered_index(conj_gp.g);

	INFO("Resolve component " << component_id);

	string graph_name = ConstructComponentName("graph_", component_id).c_str();
  // FIXME: Use path_utils
	string component_name = cfg::get().output_dir + "graph_components/"
			+ graph_name;

	ScanWithClusteredIndex(component_name, conj_gp, clustered_index);

	component_statistics(conj_gp, component_id, clustered_index);

	conj_graph_pack resolved_gp(cfg::get().K,
                              cfg::get().output_dir,
                              genome,
                              cfg::get().pos.max_single_gap, cfg::get().pos.careful_labeling);
	string sub_dir = "resolve_components/";

  // FIXME: Use path_utils
	string resolved_name = cfg::get().output_dir + "resolve_components"
			+ "/resolve_" + graph_name + "/";
	make_dir(resolved_name);

	WriteGraphPack(conj_gp, cfg::get().output_dir + sub_dir + graph_name	+ "_2_unresolved.dot");
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
	WARN("Jump resolver not alailable");

//	VERIFY(cfg::get().andrey_params.);
//	resolve_repeats_ml(gp, index, cfg::get().output_dir + "jump_resolve/",
//			cfg::get().andrey_params,
//			boost::optional<const paired_info_index&>(jump_index));
}

void prepare_jump_index(const Graph& g, const paired_info_index& raw_jump_index,
		paired_info_index& jump_index) {
	JumpingEstimator<Graph> estimator(raw_jump_index);
	paired_info_index clustered_jump_index(g);
	estimator.Estimate(clustered_jump_index);

	JumpingNormilizerFunction<Graph> nf(g, *cfg::get().ds.RL, 500);
	PairedInfoNormalizer<Graph> normalizer(nf);
	paired_info_index normalized_jump_index(g);
	normalizer.FillNormalizedIndex(clustered_jump_index, normalized_jump_index);

	JumpingPairInfoChecker<Graph> filter(g, 300, 100, 100);
	filter.Filter(normalized_jump_index, jump_index);
}

void prepare_scaffolding_index(conj_graph_pack& gp, paired_info_index& paired_index,
                                paired_info_index& clustered_index) {
    double is_var = *cfg::get().ds.is_var;
    size_t delta = size_t(is_var);
    size_t linkage_distance = size_t(cfg::get().de.linkage_distance_coeff * is_var);
    GraphDistanceFinder<Graph> dist_finder(gp.g, *cfg::get().ds.IS, *cfg::get().ds.RL, delta);

    size_t max_distance = size_t(cfg::get().de.max_distance_coeff * is_var);
    INFO("Symmetry trick");
    paired_info_index symmetric_index(gp.g);
    PairedInfoSymmetryHack<Graph> hack(gp.g, paired_index);
    hack.FillSymmetricIndex(symmetric_index);

    boost::function<double(int)> weight_function;

    INFO("Retaining insert size distribution for it");
    InsertSizeHistogramCounter<conj_graph_pack>::hist_type insert_size_hist = cfg::get().ds.hist;
    WeightDEWrapper wrapper(insert_size_hist, *cfg::get().ds.IS);
    INFO("Weight Wrapper Done");
    weight_function = boost::bind(&WeightDEWrapper::CountWeight, wrapper, _1);

    PairedInfoNormalizer<Graph>::WeightNormalizer normalizing_f;
    if (cfg::get().ds.single_cell) {
        normalizing_f = &TrivialWeightNormalization<Graph>;
    } else {
        //todo reduce number of constructor params
        PairedInfoWeightNormalizer<Graph> weight_normalizer(gp.g,
                *cfg::get().ds.IS, *cfg::get().ds.is_var, *cfg::get().ds.RL,
                gp.k_value, *cfg::get().ds.avg_coverage);
        normalizing_f = boost::bind(
                &PairedInfoWeightNormalizer<Graph>::NormalizeWeight,
                weight_normalizer, _1);
    }
    PairedInfoNormalizer<Graph> normalizer(normalizing_f);
    INFO("Normalizer Done");

    PairInfoWeightFilter<Graph> filter(gp.g, 0.);
    INFO("Weight Filter Done");

    const AbstractDistanceEstimator<Graph>& estimator =
            SmoothingDistanceEstimator<Graph>(gp.g, symmetric_index,
                    dist_finder, weight_function, linkage_distance, max_distance,
                    cfg::get().ade.threshold,
                    cfg::get().ade.range_coeff,
                    cfg::get().ade.delta_coeff, cfg::get().ade.cutoff,
                    cfg::get().ade.min_peak_points,
                    cfg::get().ade.inv_density,
                    cfg::get().ade.percentage,
                    cfg::get().ade.derivative_threshold,
                    true);
    INFO("Starting SMOOTHING distance estimator");
    estimate_with_estimator(gp.g, estimator, normalizer, filter, clustered_index);
}

void resolve_repeats() {
	Sequence genome = cfg::get().developer_mode ? cfg::get().ds.reference_genome : Sequence();

	conj_graph_pack conj_gp(cfg::get().K,
                          cfg::get().output_dir,
                          genome, cfg::get().pos.max_single_gap,
                          cfg::get().pos.careful_labeling, /*use_inner_ids*/
                          !cfg::get().developer_mode);
	paired_info_index paired_index(conj_gp.g, cfg::get().online_clust_rad);
	paired_info_index clustered_index(conj_gp.g);
	if (!cfg::get().developer_mode) {
		//Detaching edge_pos handler
		conj_gp.edge_pos.Detach();
		paired_index.Detach();
		clustered_index.Detach();
		if (!cfg::get().gap_closer_enable && !cfg::get().paired_mode) {
			conj_gp.kmer_mapper.Detach();
		}
	}

	exec_distance_estimation(conj_gp, paired_index, clustered_index);

	DEBUG("Online clusterization rad = " << cfg::get().online_clust_rad);
	if (cfg::get().developer_mode && cfg::get().pos.late_threading) {
		FillPos(conj_gp, conj_gp.genome, "10");
		FillPos(conj_gp, !conj_gp.genome, "11");
		if (!cfg::get().pos.contigs_for_threading.empty()
				&& fileExists(cfg::get().pos.contigs_for_threading)) {
			FillPos(conj_gp, cfg::get().pos.contigs_for_threading, "thr_");
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

	if (cfg::get().SAM_writer_enable && cfg::get().sw.align_before_RR) {
		SAMBeforeResolve(conj_gp);
	}

	if (!cfg::get().paired_mode
			|| cfg::get().rm == debruijn_graph::resolving_mode::rm_none) {
		OutputContigs(conj_gp.g, cfg::get().output_dir + "final_contigs.fasta");
		return;
	}

	//tSeparatedStats(conj_gp, conj_gp.genome, clustered_index);

	INFO("STAGE == Resolving Repeats");

	if (cfg::get().rm == debruijn_graph::resolving_mode::rm_split) {
		int number_of_components = 0;

//		if (cfg::get().componential_resolve) {
//			make_dir(cfg::get().output_dir + "graph_components" + "/");
//			number_of_components = PrintGraphComponents(
//					cfg::get().output_dir + "graph_components/graph_", conj_gp,
//					*cfg::get().ds.IS + 100, clustered_index);
//			INFO("number of components " << number_of_components);
//		}

		if (cfg::get().rr.symmetric_resolve) {

			PairInfoInprover<conj_graph_pack::graph_t> pi_imp(conj_gp.g);
			pi_imp.ImprovePairedInfo(clustered_index,
					cfg::get().use_multithreading, cfg::get().max_threads);
//			if (cfg::get().use_multithreading) {
//				ParalelCorrectPairedInfo(conj_gp, clustered_index, cfg::get().max_threads);
//				ParalelCorrectPairedInfo(conj_gp, clustered_index, cfg::get().max_threads);
//			} else {
//				CorrectPairedInfo(conj_gp, clustered_index, true, false);
////				CorrectPairedInfo(conj_gp, clustered_index, true, false);
////				CorrectPairedInfo(conj_gp, clustered_index, true, false);
////				CorrectPairedInfo(conj_gp, clustered_index);
////				CorrectPairedInfo(conj_gp, clustered_index);
////				CorrectPairedInfo(conj_gp, clustered_index);
////				CorrectPairedInfo(conj_gp, clustered_index);
//			}

			save_distance_filling(conj_gp, paired_index, clustered_index);

			if (cfg::get().componential_resolve) {
				make_dir(cfg::get().output_dir + "graph_components" + "/");
				number_of_components = PrintGraphComponents(
						cfg::get().output_dir + "graph_components/graph_", conj_gp,
						*cfg::get().ds.IS + 100, clustered_index);
				INFO("number of components " << number_of_components);
			}

			conj_graph_pack resolved_gp(cfg::get().K,
                                  cfg::get().output_dir,
                                  genome,
                                  cfg::get().pos.max_single_gap,
                                  cfg::get().pos.careful_labeling);
			resolved_gp.index.Detach();

			EdgeLabelHandler<conj_graph_pack::graph_t> labels_after(
					resolved_gp.g, conj_gp.g);

			process_resolve_repeats(conj_gp, clustered_index, resolved_gp,
					"graph", labels_after, "", true, cfg::get().rr.kill_loops);

	        if (cfg::get().use_scaffolder) {
                INFO("Transfering paired information");
                PairedInfoIndex<conj_graph_pack::graph_t> resolved_graph_paired_info_cl(
                      resolved_gp.g);
                ProduceResolvedPairedInfo(conj_gp, clustered_index, resolved_gp,
                      labels_after, resolved_graph_paired_info_cl);
                PairedInfoIndex<conj_graph_pack::graph_t> resolved_graph_paired_info(
                                      resolved_gp.g);

                if (cfg::get().pe_params.param_set.scaffolder_options.cluster_info) {
                    PairedInfoIndex<conj_graph_pack::graph_t> scaff_clustered(conj_gp.g);

                    //TODO: cluster here from paired_index to scaff_clustered
                    prepare_scaffolding_index(conj_gp, paired_index, scaff_clustered);

                    PairedInfoIndex<conj_graph_pack::graph_t> resolved_graph_scaff_clustered(
                                          resolved_gp.g);
                    ProduceResolvedPairedInfo(conj_gp, scaff_clustered, resolved_gp,
                          labels_after, resolved_graph_scaff_clustered);

                    DEBUG("Resolved scaffolding index size " << resolved_graph_scaff_clustered.size());

                    INFO("Scaffolding");
                    resolve_repeats_pe(cfg::get().K, resolved_gp,
                            resolved_graph_paired_info_cl,
                            resolved_graph_scaff_clustered,
                            cfg::get().output_dir,
                            "scaffolds.fasta",
                            cfg::get().pe_params);
                }
                else  {
                    ProduceResolvedPairedInfo(conj_gp, paired_index, resolved_gp,
                            labels_after, resolved_graph_paired_info);

                    INFO("Scaffolding");
                    resolve_repeats_pe(cfg::get().K, resolved_gp,
                            resolved_graph_paired_info_cl,
                            resolved_graph_paired_info,
                            cfg::get().output_dir,
                            "scaffolds.fasta",
                            cfg::get().pe_params);
                }

                save_resolved(resolved_gp, resolved_graph_paired_info, resolved_graph_paired_info_cl);
	        }


			if (cfg::get().SAM_writer_enable && cfg::get().sw.align_after_RR) {
				SAMAfterResolve(conj_gp, resolved_gp, labels_after);
			}

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
		WARN("Jump resover unavailable so far");

	}

	if (cfg::get().rm == debruijn_graph::resolving_mode::rm_path_extend) {
	    if (cfg::get().pe_params.param_set.improve_paired_info) {
	        distance_filling(conj_gp, paired_index, clustered_index);
	    }

	    if (cfg::get().pe_params.param_set.scaffolder_options.on) {
	        if (cfg::get().pe_params.param_set.scaffolder_options.cluster_info) {
                PairedInfoIndex<conj_graph_pack::graph_t> scaff_clustered(
                                      conj_gp.g);

                //TODO: cluster here

                resolve_repeats_pe(cfg::get().K, conj_gp, clustered_index, scaff_clustered,
                        cfg::get().output_dir,
                        "scaffolds.fasta",
                        cfg::get().pe_params);
            }
            else  {
                resolve_repeats_pe(cfg::get().K, conj_gp, clustered_index, paired_index,
                        cfg::get().output_dir,
                        "scaffolds.fasta",
                        cfg::get().pe_params);
            }
	    }
	    else {
	        resolve_repeats_pe(cfg::get().K, conj_gp, clustered_index,
	                cfg::get().output_dir,
	                "final_contigs.fasta",
	                cfg::get().pe_params);
	    }

	}

	if (cfg::get().rm == debruijn_graph::resolving_mode::rm_combined
			&& !cfg::get().path_set_graph) {

        INFO("Combined resolving started");

        distance_filling(conj_gp, paired_index, clustered_index);

        conj_graph_pack resolved_gp(cfg::get().K,
                                    cfg::get().output_dir,
                                    genome,
                                    cfg::get().pos.max_single_gap,
                                    cfg::get().pos.careful_labeling);

        EdgeLabelHandler<conj_graph_pack::graph_t> labels_after(
            resolved_gp.g, conj_gp.g);

        process_resolve_repeats(conj_gp, clustered_index, resolved_gp,
            "graph", labels_after, "", false, false);


        PairedInfoIndex<conj_graph_pack::graph_t> resolved_graph_paired_info_cl(
              resolved_gp.g);
        ProduceResolvedPairedInfo(conj_gp, clustered_index, resolved_gp,
              labels_after, resolved_graph_paired_info_cl);

        PairedInfoIndex<conj_graph_pack::graph_t> resolved_graph_paired_info(
              resolved_gp.g);
        if (cfg::get().pe_params.param_set.scaffolder_options.on) {
            ProduceResolvedPairedInfo(conj_gp, paired_index, resolved_gp,
                  labels_after, resolved_graph_paired_info);
        }


        if (cfg::get().pe_params.param_set.scaffolder_options.cluster_info) {
            PairedInfoIndex<conj_graph_pack::graph_t> scaff_clustered(
                                  conj_gp.g);

            //TODO: cluster here

            PairedInfoIndex<conj_graph_pack::graph_t> resolved_graph_scaff_clustered(
                                  resolved_gp.g);
            ProduceResolvedPairedInfo(conj_gp, scaff_clustered, resolved_gp,
                  labels_after, resolved_graph_scaff_clustered);

            INFO("Scaffolding");
            resolve_repeats_pe(cfg::get().K, resolved_gp,
                    resolved_graph_paired_info_cl,
                    resolved_graph_scaff_clustered,
                    cfg::get().output_dir,
                    "scaffolds.fasta",
                    cfg::get().pe_params);
        }
        else  {
            ProduceResolvedPairedInfo(conj_gp, paired_index, resolved_gp,
                    labels_after, resolved_graph_paired_info);

            INFO("Scaffolding");
            resolve_repeats_pe(cfg::get().K, resolved_gp,
                    resolved_graph_paired_info_cl,
                    resolved_graph_paired_info,
                    cfg::get().output_dir,
                    "scaffolds.fasta",
                    cfg::get().pe_params);
        }

        if (cfg::get().run_mode) {
            save_resolved(resolved_gp, resolved_graph_paired_info, resolved_graph_paired_info_cl);
        }

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

