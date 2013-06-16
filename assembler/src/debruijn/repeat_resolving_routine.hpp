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
#include "repeat_resolving.hpp"
#include "distance_estimation_routine.hpp"
#include "io/careful_filtering_reader_wrapper.hpp"
#include "io/is_corrupting_wrapper.hpp"
#include "resolved_pair_info.hpp"
#include "graph_construction.hpp"
#include "debruijn_stats.hpp"
#include "de/distance_estimation.hpp"
#include "omni/omni_utils.hpp"

#include "omni/loop_killer.hpp"
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

typedef io::CarefulFilteringReaderWrapper<io::SingleRead> CarefulFilteringStream;

namespace debruijn_graph {

void resolve_repeats(PairedReadStream& stream, const Sequence& genome);

template<class gp_t>
void WriteGraphPack(gp_t& gp, const string& file_name) {
	ofstream filestr(file_name);
	CompositeGraphColorer<typename gp_t::graph_t> colorer(
			new FixedColorer<typename gp_t::graph_t::VertexId>("white"),
			new PositionsEdgeColorer<typename gp_t::graph_t>(gp.g,
					gp.edge_pos));

	EdgeQuality<typename gp_t::graph_t> edge_qual(gp.g, gp.index,
			gp.kmer_mapper, gp.genome);
	total_labeler_graph_struct graph_struct(gp.g, &gp.int_ids, &gp.edge_pos);
	total_labeler tot_lab(&graph_struct);
	CompositeLabeler<Graph> labeler(tot_lab, edge_qual);
	DotGraphPrinter<typename gp_t::graph_t> g_print(gp.g, labeler, colorer, " ",
			filestr);
	SimpleGraphVisualizer<typename gp_t::graph_t> gv(gp.g, g_print);
	gv.Visualize();
}

void SaveResolved(conj_graph_pack& resolved_gp,
		PairedIndexT& resolved_graph_paired_info,
		PairedIndexT& resolved_graph_paired_info_cl) {

	if (cfg::get().make_saves) {
		string p = path::append_path(cfg::get().output_saves, "split_resolved");
		INFO("Saving current state to " << p);
		PrintAll(p, resolved_gp, resolved_graph_paired_info,
				resolved_graph_paired_info_cl);
		write_lib_data(p);
	}
}

template<class graph_pack>
void FindDistanceFromRepeats(graph_pack& gp,
		EdgeLabelHandler<typename graph_pack::graph_t>& labels_after,
		map<EdgeId, pair<size_t, size_t> >& distance_to_repeats_end) {
	set<EdgeId> not_unique;
	for (auto iter = gp.g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
//		VertexId start = gp.g.EdgeStart(*iter);
//		VertexId end = gp.g.EdgeEnd(*iter);
		//Graph topology implied repeats
//		if (((gp.g.CheckUniqueIncomingEdge(end) && !gp.g.IsDeadEnd(end)) ||( gp.g.CheckUniqueOutgoingEdge(start)  && !gp.g.IsDeadStart(start) )) && (gp.g.length(*iter) <  cfg::get().rr.max_repeat_length))
//			not_unique.insert(*iter);
//		//Split-based repeats
//		else
		if (labels_after.edge_inclusions.find(*iter)
				!= labels_after.edge_inclusions.end()
				&& labels_after.edge_inclusions[*iter].size() > 1)
			not_unique.insert(*iter);
	}
	set<set<EdgeId> > components;
	while (!not_unique.empty()) {
		set<EdgeId> wfs_set;
		wfs_set.insert(*not_unique.begin());
		set<EdgeId> component;
		while (!wfs_set.empty()) {
			for (auto iter = wfs_set.begin(); iter != wfs_set.end();) {
				component.insert(*iter);
				VertexId start = gp.g.EdgeStart(*iter);
				VertexId end = gp.g.EdgeEnd(*iter);
				vector<EdgeId> next = gp.g.IncomingEdges(start);
				for (auto e_iter = next.begin(); e_iter != next.end(); e_iter++)
					if (not_unique.find(*e_iter) != not_unique.end()
							&& component.find(*e_iter) == component.end())
						wfs_set.insert(*e_iter);
				next = gp.g.OutgoingEdges(end);
				for (auto e_iter = next.begin(); e_iter != next.end(); e_iter++)
					if (not_unique.find(*e_iter) != not_unique.end()
							&& component.find(*e_iter) == component.end())
						wfs_set.insert(*e_iter);
				set<EdgeId>::iterator new_iter = iter;
				new_iter++;
				wfs_set.erase(iter);
				iter = new_iter;

			}
		}
		DEBUG("not_unique_component");
		for (auto iter = component.begin(); iter != component.end(); ++iter) {
			DEBUG(gp.g.int_id(*iter));
			not_unique.erase(*iter);
		}
		components.insert(component);
		FillComponentDistances(component, distance_to_repeats_end, gp);
		component.clear();
	}
	string p = path::append_path(cfg::get().output_saves,
			"distance_from_repeats");
	if (cfg::get().make_saves)
		SaveComponents(p, components, gp);
}

template<class graph_pack>
void FillComponentDistances(set<EdgeId>& component,
		map<EdgeId, pair<size_t, size_t> > & distances_map, graph_pack& gp) {
	map<EdgeId, pair<size_t, size_t> > component_map;
	map<VertexId, pair<int, int> > vertex_map;
	//longest pathes from start and end of edge resp to first base(backward and forward resp) not in repeat;
	for (auto iter = component.begin(); iter != component.end(); iter++) {
		component_map.insert(
				std::make_pair(*iter, std::make_pair(1000000, 1000000)));
		VertexId start = gp.g.EdgeStart(*iter);
		VertexId end = gp.g.EdgeEnd(*iter);
		vertex_map.insert(
				std::make_pair(start, std::make_pair(-1000000, -1000000)));
		vertex_map.insert(
				std::make_pair(end, std::make_pair(-1000000, -1000000)));
	}
	for (auto iter = component.begin(); iter != component.end(); iter++) {
		VertexId start = gp.g.EdgeStart(*iter);
		VertexId end = gp.g.EdgeEnd(*iter);
		vector<EdgeId> next = gp.g.IncomingEdges(start);
		for (auto e_iter = next.begin(); e_iter != next.end(); e_iter++)
			if (component.find(*e_iter) == component.end())
				vertex_map[start].first = 0;
		if (next.size() == 0)
			vertex_map[start].first = 0;

		next = gp.g.OutgoingEdges(end);
		for (auto e_iter = next.begin(); e_iter != next.end(); e_iter++)
			if (component.find(*e_iter) == component.end())
				vertex_map[end].second = 0;
		if (next.size() == 0)
			vertex_map[end].second = 0;
	}
	pair<set<EdgeId>, set<EdgeId> > used;
	for (size_t j = 0; j < 10000; j++) {
		bool changed = false;
		for (auto iter = component.begin(); iter != component.end(); iter++) {
			VertexId start = gp.g.EdgeStart(*iter);
			VertexId end = gp.g.EdgeEnd(*iter);
			int len = gp.g.length(*iter);
			if (used.first.find(*iter) == used.first.end()) {
				if (vertex_map[start].first >= 0
						&& vertex_map[start].first + len
								> vertex_map[end].first) {
					vertex_map[end].first = vertex_map[start].first + len;
					changed = true;
					used.first.insert(*iter);
				}
			}

			if (used.second.find(*iter) == used.second.end()) {
				if (vertex_map[end].second >= 0
						&& vertex_map[end].second + len
								> vertex_map[start].second) {
					vertex_map[start].second = vertex_map[end].second + len;
					changed = true;
					used.second.insert(*iter);
				}
			}
		}
		if (!changed)
			break;
	}
	for (auto iter = component.begin(); iter != component.end(); iter++) {
		VertexId start = gp.g.EdgeStart(*iter);
		VertexId end = gp.g.EdgeEnd(*iter);
		component_map[*iter] = std::make_pair(vertex_map[start].first,
				vertex_map[end].second);
	}
//	for(size_t j = 0; j < 1000; j++){
//		bool changed = false;
//		for (auto iter = component.begin(); iter != component.end(); iter ++) {
//			size_t max_incoming_path = -1000;
//			size_t max_outgoing_path = -1000;
//			VertexId start = gp.g.EdgeStart(*iter);
//			VertexId end = gp.g.EdgeEnd(*iter);
//			vector<EdgeId> next = gp.g.IncomingEdges(start);
//			for (auto e_iter = next.begin(); e_iter != next.end(); e_iter++)
//				if (component_map.find(*e_iter) != component_map.end())
//					max_incoming_path = max(max_incoming_path, component_map[*e_iter].first + gp.g.length(*e_iter));
//			next = gp.g.OutgoingEdges(end);
//			for (auto e_iter = next.begin(); e_iter != next.end(); e_iter++)
//				if (component_map.find(*e_iter) != component_map.end())
//					max_outgoing_path = min(max_outgoing_path, component_map[*e_iter].second + gp.g.length(*e_iter));
//			if (max_incoming_path < component_map[*iter].first){
//				changed = true;
//				component_map[*iter].first = max_incoming_path;
//			}
//			if (max_outgoing_path < component_map[*iter].second){
//				changed = true;
//				component_map[*iter].second = max_outgoing_path;
//			}
//		}
//		if (!changed)
//			break;
//	}
	DEBUG("not_unique_component");
	for (auto iter = component.begin(); iter != component.end(); ++iter) {
		DEBUG(
				gp.g.int_id(*iter)<<" len "<< gp.g.length(*iter) <<" :  "<<component_map[*iter].first <<" " << component_map[*iter].second);
		distances_map[*iter] = component_map[*iter];
	}
}

//TODO Move to graphio if saves needed;
template<class graph_pack>
void SaveComponents(string file_name, set<set<EdgeId> >& components,
		graph_pack& gp) {
	FILE* file = fopen((file_name + ".rep").c_str(), "w");
	fprintf(file, "%zu\n", components.size());
	for (auto iter = components.begin(); iter != components.end(); ++iter) {
		fprintf(file, "%zu\n", iter->size());
		for (auto j_iter = iter->begin(); j_iter != iter->end(); ++j_iter) {
			fprintf(file, "%zu ", gp.int_ids.ReturnIntId(*j_iter));
		}
		fprintf(file, "\n");
	}
	fclose(file);
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
		const PairedInfoIndexT<typename graph_pack::graph_t>& clustered_index,
		graph_pack& resolved_gp,
		EdgeLabelHandler<typename graph_pack::graph_t>& labels_after,
		PairedInfoIndexT<typename graph_pack::graph_t>& resolved_graph_paired_info) {
	INFO("Generating paired info for resolved graph");
	ResolvedGraphPairInfoCounter<typename graph_pack::graph_t> resolved_graph_paired_info_counter(
			origin_gp.g, clustered_index, resolved_gp.g, labels_after);
	resolved_graph_paired_info_counter.FillResolvedGraphPairedInfo(
			resolved_graph_paired_info);
	DEBUG("Generating paired info for resolved graph complete");
}

template<class graph_pack>
void SaveResolvedPairedInfo(graph_pack& resolved_gp,
		PairedInfoIndexT<typename graph_pack::graph_t> resolved_graph_paired_info,
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
		PrintWithClusteredIndex(rr_filename, resolved_gp,
				resolved_graph_paired_info);
		DEBUG("Saved");
	}
}

bool try_load_paired_index(conj_graph_pack& gp, PairedIndexT& index,
		const string& name, path::files_t* used_files = 0) {

	string p = path::append_path(cfg::get().load_from, name);

	FILE* file = fopen((p + ".prd").c_str(), "r");
	if (file == NULL) {
		return false;
	}
	fclose(file);

	if (used_files != 0)
		used_files->push_back(p);

	index.Clear();
	ScannerTraits<conj_graph_pack::graph_t>::Scanner scanner(gp.g, gp.int_ids);
	ScanPairedIndex(p, scanner, index);

	return true;

}

template<class graph_pack>
void RemapMaskedMismatches(graph_pack& resolved_gp, graph_pack& origin_gp,
        const io::SequencingLibrary<debruijn_config::DataSetData> &lib,
		EdgeLabelHandler<typename graph_pack::graph_t>& labels_after,
		map<EdgeId, pair<size_t, size_t> >& distance_to_repeats_end) {
	size_t Ncount = 0;
	//FILE* file = fopen(("multipicities.tmp"), "w");
	int not_masked = 0;
	for (auto iter = origin_gp.g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
//		size_t len = origin_gp.g.length(*iter) + origin_gp.g.k();
		size_t multiplicity = labels_after.edge_inclusions[*iter].size();
//		fprintf(file, "%d %d \n",origin_gp.int_ids.ReturnIntId(*iter), multiplicity);
		if (multiplicity > 0) {
//			INFO(origin_gp.g.int_id(*iter));
//			INFO(origin_gp.g.int_id(origin_gp.g.conjugate(*iter)));
			const vector<
					typename MismatchMasker<typename graph_pack::graph_t>::MismatchInfo> mismatches =
					origin_gp.mismatch_masker.mismatch_map[*iter];
			DEBUG(mismatches.size());
			const vector<
					typename MismatchMasker<typename graph_pack::graph_t>::MismatchInfo> rc_mismatches =
					origin_gp.mismatch_masker.mismatch_map[origin_gp.g.conjugate(
							*iter)];
			DEBUG(rc_mismatches.size());
			if (mismatches.size() != rc_mismatches.size()) {
				WARN(mismatches.size() <<" /// " << rc_mismatches.size());
			}

			for (size_t i = 0; i < mismatches.size(); i++) {
//TODO:: cutoff selection!
				vector<pair<EdgeId, size_t> > resolved_positions =
						labels_after.resolvedPositions(*iter,
								mismatches[i].position);
				double cutoff = 0.5;
				if ((origin_gp.g.length(*iter) > size_t(lib.data().mean_insert_size)
						&& multiplicity > 1)
						|| (distance_to_repeats_end[*iter].first
								+ mismatches[i].position > size_t(lib.data().mean_insert_size)
								&& distance_to_repeats_end[*iter].second
										+ origin_gp.g.length(*iter)
										- mismatches[i].position
										> size_t(lib.data().mean_insert_size)))
					cutoff /= (4); // /cfg::get().mismatch_ratio);
				else {
					cutoff *= 1.5; //* cfg::get().mismatch_ratio;
				}
				bool cfg_corrector_on = true;
				if (cfg_corrector_on
						&& !(distance_to_repeats_end[*iter].first
								+ mismatches[i].position > size_t(lib.data().mean_insert_size)
								&& distance_to_repeats_end[*iter].second
										+ origin_gp.g.length(*iter)
										- mismatches[i].position
										> size_t(lib.data().mean_insert_size))) {
					not_masked++;
					origin_gp.mismatch_masker.mismatch_map[*iter][i].ratio = 0;
					continue;
				}
				map<EdgeId, int> diff_res;
				for (auto it = resolved_positions.begin();
						it < resolved_positions.end(); it++)
					if (diff_res.find(it->first) == diff_res.end())
						diff_res.insert(make_pair(it->first, 1));
					else
						diff_res[it->first]++;
				int real_count = 0;
				for (size_t j = 0; j < resolved_positions.size(); j++) {
					double real_multiplicity = origin_gp.g.coverage(*iter)
							/ resolved_gp.g.coverage(
									resolved_positions[j].first);
					if (real_multiplicity
							* diff_res[resolved_positions[j].first]
							* mismatches[i].ratio > cutoff)
						real_count++;
				}
				for (size_t j = 0; j < resolved_positions.size(); j++) {
//					if (origin_gp.g.int_id(*iter) >= 10006349 && origin_gp.g.int_id(*iter) <= 10006352) {
//						INFO(origin_gp.g.int_id(*iter) << " cutoff: " << cutoff);
//						INFO(origin_gp.g.coverage(*iter) <<" // " << resolved_gp.g.coverage(resolved_positions[j].first));
//						INFO("loop coefficient" << diff_res[resolved_positions[j].first]);
//					}
					double real_multiplicity = origin_gp.g.coverage(*iter)
							/ resolved_gp.g.coverage(
									resolved_positions[j].first);

					if (real_multiplicity
							* diff_res[resolved_positions[j].first]
							* mismatches[i].ratio > cutoff && real_count <= 5) {
						DEBUG(
								origin_gp.g.int_id(*iter) << " position after: "<< resolved_positions[j].second << "position before: "<< mismatches[i].position << " edge length:  "<< origin_gp.g.length(*iter));
						DEBUG(
								distance_to_repeats_end[*iter].first<< " " << distance_to_repeats_end[*iter].second);

						resolved_gp.mismatch_masker.insert(
								resolved_positions[j].first,
								resolved_positions[j].second,
								real_multiplicity * mismatches[i].ratio,
								mismatches[i].counts, cutoff);
						Ncount++;
					}
				}

				if (cutoff > mismatches[i].ratio)
					origin_gp.mismatch_masker.mismatch_map[*iter][i].ratio = 0;
			}
		}
	}

	INFO("masked "<< Ncount << " potential mismatches masked");
	INFO(" "<< not_masked<< " potential mismatches left for corrector");
}

template<class graph_pack>
void process_resolve_repeats(graph_pack& origin_gp,
        const io::SequencingLibrary<debruijn_config::DataSetData> &lib,
		const PairedInfoIndexT<typename graph_pack::graph_t>& clustered_index,
		graph_pack& resolved_gp, const string& graph_name,
		EdgeLabelHandler<typename graph_pack::graph_t>& labels_after,
		const string& subfolder = "", bool output_contigs = true,
		bool kill_loops = true) {
	typename graph_pack::graph_t &g = origin_gp.g;
	const PairedInfoIndexT<typename graph_pack::graph_t> &pii = clustered_index;
	if (cfg::get().cut_bad_connections)
		BadConnectionCutter<typename graph_pack::graph_t>(g, pii).CutConnections();
//	EdgeLabelHandler<typename graph_pack::graph_t> labels_after(resolved_gp.g,
//			origin_gp.g);
//	ProduceLongEdgesStat( origin_gp,  clustered_index);
	if (cfg::get().compute_paths_number)
		GenerateMatePairStats(origin_gp, clustered_index, lib.data().insert_size_deviation);
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

//    CleanIsolated(origin_gp);
	ResolveRepeats(origin_gp.g, lib, origin_gp.int_ids, clustered_index,
			origin_gp.edge_pos, resolved_gp.g, resolved_gp.int_ids,
			resolved_gp.edge_pos,
			cfg::get().output_dir + subfolder + "resolve_" + graph_name + "/",
			labels_after, cfg::get().developer_mode);

	//Generating paired info for resolved graph
//    PairedInfoIndexT<typename graph_pack::graph_t> resolved_graph_paired_info(resolved_gp.g);
//		ProduceResolvedPairedInfo(origin_gp, clustered_index, resolved_gp, labels_after, resolved_graph_paired_info);
//		SaveResolvedPairedInfo(resolved_gp, resolved_graph_paired_info, graph_name + "_resolved", subfolder);
	//Paired info for resolved graph generated

	if (cfg::get().output_nonfinal_contigs && output_contigs) {
		OutputContigs(resolved_gp.g,
				cfg::get().output_dir + "after_rr_before_simplify" + postfix);
		OutputContigs(origin_gp.g,
				cfg::get().output_dir + "before_resolve.fasta");
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
		PairedInfoIndexT<typename graph_pack::graph_t> resolved_cleared_graph_paired_info_before(
				resolved_gp.g);
		ProduceResolvedPairedInfo(origin_gp, clustered_index, resolved_gp,
				labels_after, resolved_cleared_graph_paired_info_before);
	}

	INFO("SUBSTAGE == Clearing resolved graph");

	omnigraph::Compressor<typename graph_pack::graph_t> compressor(
			resolved_gp.g);
	compressor.CompressAllVertices();

//	    omnigraph::StrGraphLabeler<typename graph_pack::graph_t> str_labeler(resolved_gp.g);
//	omnigraph::WriteSimple(resolved_gp.g, str_labeler,
//				cfg::get().output_dir + subfolder + graph_name + "_resolved.dot",
//				"no_repeat_graph");

	size_t iters = 3; // TODO Constant 3? Shouldn't it be taken from config?
	map<EdgeId, pair<size_t, size_t> > distance_to_repeats_end;
	FindDistanceFromRepeats(origin_gp, labels_after, distance_to_repeats_end);
	RemapMaskedMismatches(resolved_gp, origin_gp, lib, labels_after,
			distance_to_repeats_end);
	for (size_t i = 0; i < iters; ++i) {
		INFO(
				"Tip clipping iteration " << i << " (0-indexed) out of " << iters << ":");
//		ClipTipsForResolver(resolved_gp.g);
		ClipTips(resolved_gp.g, cfg::get().simp.tc, lib.data().read_length, /*max_coverage*/
				0.);

		//PairedInfoIndexT<typename graph_pack::graph_t> resolved_cleared_graph_paired_info_before(
		//resolved_gp.g);
		//ProduceResolvedPairedInfo(origin_gp, clustered_index, resolved_gp,
		//labels_after, resolved_cleared_graph_paired_info_before);

//		INFO("Erroneous remove "<<i);
//        BulgeRemoveWrap      (resolved_gp.g);
//		FinalRemoveErroneousEdges(resolved_gp.g, edge_remover);

		//PairedInfoIndexT<typename graph_pack::graph_t> resolved_cleared_graph_paired_info_before(
		//resolved_gp.g);

		//ProduceResolvedPairedInfo(origin_gp, clustered_index, resolved_gp,
		//labels_after, resolved_cleared_graph_paired_info_before);
//        RemoveRelativelyLowCoverageEdges(resolved_gp.g);
//		omnigraph::WriteSimple(resolved_gp.g, tot_labeler_after,
//
//		cfg::get().output_dir + subfolder + ToString(i) + "b_4_cleared.dot",
//				"no_repeat_graph");
	}


    if (kill_loops) {
        SimpleLoopKiller<typename graph_pack::graph_t> lk(resolved_gp.g,
                cfg::get().rr.max_repeat_length, 6);
        lk.KillAllLoops();
    }

    OutputMaskedContigs(origin_gp.g,
			cfg::get().output_dir + "before_resolve_masked.fasta",
			origin_gp.mismatch_masker);

	DEBUG("Clearing resolved graph complete");


	//Generating paired info for resolved graph
	PairedInfoIndexT<typename graph_pack::graph_t> resolved_cleared_graph_paired_info(
			resolved_gp.g);

	ProduceResolvedPairedInfo(origin_gp, clustered_index, resolved_gp,
			labels_after, resolved_cleared_graph_paired_info);
	SaveResolvedPairedInfo(resolved_gp, resolved_cleared_graph_paired_info,
			graph_name + "_resolved_cleared", subfolder);
	//Paired info for resolved graph generated

	DEBUG("Output Contigs");

	if (output_contigs) {
		OutputMaskedContigs(resolved_gp.g,
				cfg::get().output_dir + "final_contigs_masked.fasta",
				resolved_gp.mismatch_masker);
		cfg::get_writable().final_contigs_file = cfg::get().output_dir
				+ "final_contigs.fasta";
		OutputContigs(resolved_gp.g,
				cfg::get().output_dir + "final_contigs_unmasked.fasta");
		if (cfg::get().paired_mode) {
			INFO("Outputting final masked contigs: ");
			OutputMaskedContigs(resolved_gp.g,
					cfg::get().output_dir + "final_contigs.fasta",
					resolved_gp.mismatch_masker, false, 0,
					cfg::get().cut_bad_connections);
		} else {
			INFO("Outputting final NONmasked contigs: ");
			OutputContigs(resolved_gp.g,
					cfg::get().output_dir + "final_contigs.fasta", 0,
					cfg::get().cut_bad_connections);
		}
	}
	OutputCutContigs(resolved_gp.g, cfg::get().output_dir + "cut.fasta");

	if (cfg::get().output_pictures) {
		omnigraph::WriteSimple(resolved_gp.g, tot_labeler_after,
				cfg::get().output_dir + subfolder + graph_name
						+ "_4_cleared.dot", "no_repeat_graph");
		string file_str = cfg::get().output_dir + subfolder + graph_name
				+ "_4_cleared_colored.dot";
		ofstream filestr(file_str.c_str());
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
        const io::SequencingLibrary<debruijn_config::DataSetData> &lib,
		PairedInfoIndexT<typename graph_pack::graph_t>& clustered_index,
		graph_pack& resolved_gp, const string& graph_name,
		const string& subfolder = "", bool output_contigs = true) {

	EdgeLabelHandler<typename graph_pack::graph_t> labels_after(resolved_gp.g,
			origin_gp.g);

	process_resolve_repeats(origin_gp, lib, clustered_index, resolved_gp, graph_name,
			labels_after, subfolder, output_contigs);
}

template<class graph_pack>
void component_statistics(graph_pack & conj_gp,
        const io::SequencingLibrary<debruijn_config::DataSetData> &lib,
        int component_id,
		PairedInfoIndexT<typename graph_pack::graph_t>& clustered_index) {

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
		if (conj_gp.g.length(*iter) > size_t(lib.data().mean_insert_size) + 100) {
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

void resolve_conjugate_component(int component_id, const Sequence& genome, const io::SequencingLibrary<debruijn_config::DataSetData> &lib) {
	conj_graph_pack conj_gp(cfg::get().K, cfg::get().output_dir, genome,
			cfg::get().pos.max_single_gap, cfg::get().pos.careful_labeling);
	PairedIndexT paired_index(conj_gp.g/*, 5.*/);
	PairedIndexT clustered_index(conj_gp.g);

	INFO("Resolve component " << component_id);

	string graph_name = ConstructComponentName("graph_", component_id).c_str();
	// FIXME: Use path_utils
	string component_name = cfg::get().output_dir + "graph_components/"
			+ graph_name;

	ScanWithClusteredIndex(component_name, conj_gp, clustered_index);

	component_statistics(conj_gp, lib, component_id, clustered_index);

	conj_graph_pack resolved_gp(cfg::get().K, cfg::get().output_dir, genome,
			cfg::get().pos.max_single_gap, cfg::get().pos.careful_labeling);
	string sub_dir = "resolve_components/";

	// FIXME: Use path_utils
	string resolved_name = cfg::get().output_dir + "resolve_components"
			+ "/resolve_" + graph_name + "/";
	make_dir(resolved_name);

	WriteGraphPack(conj_gp,
			cfg::get().output_dir + sub_dir + graph_name + "_2_unresolved.dot");
	process_resolve_repeats(conj_gp, lib, clustered_index, resolved_gp, graph_name,
			sub_dir, false);
}

void resolve_nonconjugate_component(int /*component_id*/, const Sequence& /*genome*/) {
//	nonconj_graph_pack nonconj_gp(genome);
//  PairedInfoIndexT<nonconj_graph_pack::graph_t> clustered_index(nonconj_gp.g);
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

void resolve_with_jumps(conj_graph_pack& /*gp*/, PairedInfoIndexT<Graph>& /*index*/,
		const PairedIndexT& /*jump_index*/) {
	WARN("Jump resolver not alailable");

//	VERIFY(cfg::get().andrey_params.);
//	resolve_repeats_ml(gp, index, cfg::get().output_dir + "jump_resolve/",
//			cfg::get().andrey_params,
//      boost::optional<const PairedIndexT>(jump_index));
}

void prepare_jump_index(const Graph& g, const PairedIndexT& raw_jump_index,
		PairedIndexT& jump_index) {
	JumpingEstimator<Graph> estimator(raw_jump_index);
	PairedIndexT clustered_jump_index(g);
	estimator.Estimate(clustered_jump_index);

	JumpingNormalizerFunction<Graph> nf(g, cfg::get().ds.RL(), 500);
	PairedInfoNormalizer<Graph> normalizer(nf);
	PairedIndexT normalized_jump_index(g);
	normalizer.FillNormalizedIndex(clustered_jump_index, normalized_jump_index);

	JumpingPairInfoChecker<Graph> filter(g, 300, 100, 100);
	filter.Filter(normalized_jump_index, jump_index);
}

void prepare_scaffolding_index(conj_graph_pack& gp,
        const io::SequencingLibrary<debruijn_config::DataSetData> &lib,
        PairedIndexT& paired_index,
		PairedIndexT& clustered_index) {
	double is_var = lib.data().insert_size_deviation;
	size_t delta = size_t(is_var);
	size_t linkage_distance = size_t(
			cfg::get().de.linkage_distance_coeff * is_var);
	GraphDistanceFinder<Graph> dist_finder(gp.g, (size_t)math::round(lib.data().mean_insert_size),
	        lib.data().read_length, delta);
	size_t max_distance = size_t(cfg::get().de.max_distance_coeff * is_var);
	boost::function<double(int)> weight_function;
	INFO("Retaining insert size distribution for it");
	map<int, size_t> insert_size_hist = cfg::get().ds.hist();
	WeightDEWrapper wrapper(insert_size_hist, size_t(lib.data().mean_insert_size));
	INFO("Weight Wrapper Done");
	weight_function = boost::bind(&WeightDEWrapper::CountWeight, wrapper, _1);

	PairedInfoNormalizer<Graph>::WeightNormalizer normalizing_f;
	if (cfg::get().ds.single_cell) {
		normalizing_f = &TrivialWeightNormalization<Graph>;
	} else {
		//todo reduce number of constructor params
	    //TODO: apply new system
		PairedInfoWeightNormalizer<Graph> weight_normalizer(gp.g,
                (size_t)math::round(lib.data().mean_insert_size), lib.data().insert_size_deviation, lib.data().read_length,
				gp.k_value, lib.data().average_coverage);
		normalizing_f = boost::bind(
				&PairedInfoWeightNormalizer<Graph>::NormalizeWeight,
				weight_normalizer, _1, _2, _3);
	}
	PairedInfoNormalizer<Graph> normalizer(normalizing_f);
	INFO("Normalizer Done");

	PairInfoWeightFilter<Graph> filter(gp.g, 0.);
	INFO("Weight Filter Done");

	const AbstractDistanceEstimator<Graph>& estimator =
			SmoothingDistanceEstimator<Graph>(gp.g, paired_index, dist_finder,
					weight_function, linkage_distance, max_distance,
					cfg::get().ade.threshold, cfg::get().ade.range_coeff,
					cfg::get().ade.delta_coeff, cfg::get().ade.cutoff,
					cfg::get().ade.min_peak_points, cfg::get().ade.inv_density,
					cfg::get().ade.percentage,
					cfg::get().ade.derivative_threshold, true);
	estimate_with_estimator(gp.g, estimator, normalizer, filter,
			clustered_index);
}

void resolve_repeats_by_coverage(conj_graph_pack& conj_gp, size_t insert_size, std::vector< PathInfo<Graph> >& filteredPaths, 
				PairedIndexT& clustered_index,
				EdgeQuality<Graph>& quality_labeler ) {


	DeBruijnEdgeIndex<conj_graph_pack::graph_t, runtime_k::RtSeq> kmerIndex(conj_gp.index.inner_index().K(), conj_gp.g, cfg::get().output_dir);
	if (cfg::get().developer_mode) {

		std::string path;
		if (cfg::get().entry_point < ws_repeats_resolving) 
			path = cfg::get().output_dir + "/saves/debruijn_kmer_index_after_construction";
		else
			path = cfg::get().load_from + "/debruijn_kmer_index_after_construction";
		bool val = LoadEdgeIndex(path, kmerIndex);
		VERIFY_MSG(val, "can not open file "+path+".kmidx");
		INFO("Updating index from graph started");
		DeBruijnEdgeIndexBuilder<runtime_k::RtSeq>().UpdateIndexFromGraph<conj_graph_pack::graph_t>(kmerIndex, conj_gp.g);
	}

	auto index = FlankingCoverage<conj_graph_pack::graph_t>(conj_gp.g, kmerIndex, 50, cfg::get().K + 1);
	EdgeLabelHandler<conj_graph_pack::graph_t> labels_after(conj_gp.g, conj_gp.g);
	auto cov_rr = CoverageBasedResolution<conj_graph_pack> (&conj_gp, cfg::get().coverage_threshold_one_list, cfg::get().coverage_threshold_match, 
			cfg::get().coverage_threshold_global, cfg::get().tandem_ratio_lower_threshold, cfg::get().tandem_ratio_upper_threshold, cfg::get().repeat_length_upper_threshold);
	cov_rr.resolve_repeats_by_coverage(index, insert_size, labels_after, quality_labeler, clustered_index, filteredPaths);

	INFO("Repeats are resolved by coverage");
}

int get_first_pe_lib_index() {
	for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
		if (cfg::get().ds.reads[i].type() == io::LibraryType::PairedEnd) {
			return i;
		}
	}
	return -1;
}

void prepare_all_scaf_libs(conj_graph_pack& conj_gp,
		vector<PairedIndexT*>& scaff_indexs, vector<size_t>& indexes) {

	vector<PairedIndexT*> cl_scaff_indexs;
	for (size_t i = 0; i < scaff_indexs.size(); ++i) {
		PairedIndexT* pe = new PairedIndexT(conj_gp.g);
		cl_scaff_indexs.push_back(pe);
		prepare_scaffolding_index(conj_gp, cfg::get().ds.reads[indexes[i]], *scaff_indexs[i], *cl_scaff_indexs[i]);
	}
	scaff_indexs.clear();
	scaff_indexs.insert(scaff_indexs.end(), cl_scaff_indexs.begin(),
			cl_scaff_indexs.end());
}

void delete_index(vector<PairedIndexT*>& index){
	for (size_t i = 0; i < index.size(); ++i){
		delete index[i];
	}
}

//Use only one pe library
void split_resolving(conj_graph_pack& conj_gp, PairedIndicesT& paired_indices,
		PairedIndicesT& clustered_indices, Sequence& genome,
		size_t pe_lib_index) {

	PairedIndexT& paired_index = paired_indices[pe_lib_index];
	PairedIndexT& clustered_index = clustered_indices[pe_lib_index];
	const io::SequencingLibrary<debruijn_config::DataSetData> &lib = cfg::get().ds.reads[pe_lib_index];
	int number_of_components = 0;
	//TODO: do we have non symmetric_resolve? can we delete this if?
	if (cfg::get().rr.symmetric_resolve) {
		if (cfg::get().componential_resolve) {
			make_dir(cfg::get().output_dir + "graph_components" + "/");
			number_of_components = PrintGraphComponents(
					cfg::get().output_dir + "graph_components/graph_", conj_gp,
					size_t(lib.data().mean_insert_size) + 100, clustered_index);
			INFO("number of components " << number_of_components);
		}

		conj_graph_pack resolved_gp(cfg::get().K, cfg::get().output_dir, genome,
				cfg::get().pos.max_single_gap, cfg::get().pos.careful_labeling);
		resolved_gp.index.Detach();

		EdgeLabelHandler<conj_graph_pack::graph_t> labels_after(resolved_gp.g,
				conj_gp.g);

		process_resolve_repeats(conj_gp, lib, clustered_index, resolved_gp, "graph",
				labels_after, "", true, cfg::get().rr.kill_loops);

		if (cfg::get().use_scaffolder) {
			vector<size_t> indexs;
			indexs.push_back(pe_lib_index);
			INFO("Transfering paired information");
			PairedInfoIndexT<conj_graph_pack::graph_t> resolved_graph_paired_info_cl(
					resolved_gp.g);
			ProduceResolvedPairedInfo(conj_gp, clustered_index, resolved_gp,
					labels_after, resolved_graph_paired_info_cl);
			PairedInfoIndexT<conj_graph_pack::graph_t> resolved_graph_paired_info(
					resolved_gp.g);
			ProduceResolvedPairedInfo(conj_gp, paired_index, resolved_gp,
					labels_after, resolved_graph_paired_info);
			vector<PairedIndexT*> pe_indexs;
			vector<PairedIndexT*> pe_scaf_indexs;
			pe_indexs.push_back(&resolved_graph_paired_info_cl);
			PairedInfoIndexT<conj_graph_pack::graph_t>* resolved_graph_scaff_clustered =
									new PairedInfoIndexT<conj_graph_pack::graph_t>(resolved_gp.g);
			if (cfg::get().pe_params.param_set.scaffolder_options.cluster_info) {
				PairedInfoIndexT<conj_graph_pack::graph_t> scaff_clustered(
						conj_gp.g);
				prepare_scaffolding_index(conj_gp, lib, paired_index,
						scaff_clustered);
				ProduceResolvedPairedInfo(conj_gp, scaff_clustered, resolved_gp,
						labels_after, *resolved_graph_scaff_clustered);
				DEBUG("Resolved scaffolding index size " << resolved_graph_scaff_clustered->size());
				pe_scaf_indexs.push_back(resolved_graph_scaff_clustered);
			} else {
				pe_scaf_indexs.push_back(&resolved_graph_paired_info);
			}
			INFO("Scaffolding");
			path_extend::resolve_repeats_pe(resolved_gp, pe_indexs,
					pe_scaf_indexs, indexs,  vector<PathInfo<Graph > >(), cfg::get().output_dir, "scaffolds.fasta", false, boost::none, false);
			SaveResolved(resolved_gp, resolved_graph_paired_info,
					resolved_graph_paired_info_cl);
			delete resolved_graph_scaff_clustered;
		}

		if (cfg::get().componential_resolve) {
			make_dir(cfg::get().output_dir + "resolve_components" + "/");
			for (int i = 0; i < number_of_components; i++) {
				resolve_conjugate_component(i + 1, genome, lib);
			}
		}
	}
}

void pe_resolving(conj_graph_pack& conj_gp, PairedIndicesT& paired_indices,	PairedIndicesT& clustered_indices, EdgeQuality<Graph>& quality_labeler) {

	vector<PairedIndexT*> pe_indexs;
	vector<PairedIndexT*> pe_scaf_indexs;
	vector<size_t> indexes;

	for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
		if (cfg::get().ds.reads[i].type() == io::LibraryType::PairedEnd
				|| cfg::get().ds.reads[i].type() == io::LibraryType::MatePairs) {
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
	OutputContigs(conj_gp.g, cfg::get().output_dir + "before_resolve.fasta");
	if (cfg::get().coverage_based_rr == true){
		int pe_lib_index = get_first_pe_lib_index();
		const io::SequencingLibrary<debruijn_config::DataSetData> &lib = cfg::get().ds.reads[pe_lib_index];
		resolve_repeats_by_coverage(conj_gp, lib.data().mean_insert_size, filteredPaths, clustered_indices[0], quality_labeler);
	}



    //LongReadStorage<Graph> long_read(conj_gp.g);
	if (cfg::get().pacbio_test_on == true){
		INFO("creating  multiindex with k = " << cfg::get().pb.pacbio_k);
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
	}
	else {
		pe_scaf_indexs.clear();
		//path_extend::resolve_repeats_pe(conj_gp, pe_indexs, pe_scaf_indexs, indexes, long_read.GetAllPaths(), cfg::get().output_dir, "final_contigs.fasta", false, boost::none);
		path_extend::resolve_repeats_pe(conj_gp, pe_indexs, pe_scaf_indexs, indexes, filteredPaths, cfg::get().output_dir, "final_contigs.fasta", false, boost::none);
	}
}

void resolve_repeats() {

	Sequence genome =
			cfg::get().developer_mode ?
					cfg::get().ds.reference_genome : Sequence();

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

	//Repeat resolving begins
	int pe_lib_index = get_first_pe_lib_index();
	INFO("STAGE == Resolving Repeats");
	if (cfg::get().ds.reads.lib_count() == 1 && pe_lib_index >= 0
			&& cfg::get().rm == debruijn_graph::resolving_mode::rm_split) {
		INFO("Split repeat resolving");
		split_resolving(conj_gp, paired_indices, clustered_indices, genome,
				pe_lib_index);
	}
	else if (cfg::get().ds.reads.lib_count() > 1 || pe_lib_index == -1
			|| cfg::get().rm
					== debruijn_graph::resolving_mode::rm_path_extend) {
		INFO("Path-Extend repeat resolving");
		pe_resolving(conj_gp, paired_indices, clustered_indices, quality_labeler);
	}
	else if (cfg::get().rm == debruijn_graph::resolving_mode::rm_rectangles) {
		INFO("Ready to run rectangles repeat resolution module");
	}
	else {
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

