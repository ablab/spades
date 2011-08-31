/*
 * launch.hpp
 *
 *  Created on: May 6, 2011
 *      Author: sergey
 */

#ifndef LAUNCH_HPP_
#define LAUNCH_HPP_

#include "io/converting_reader_wrapper.hpp"
#include "omni/visualization_utils.hpp"

//#include "debruijn_graph.hpp"
#include "omni/edge_labels_handler.hpp"
#include "omni/paired_info.hpp"
#include "graph_construction.hpp"
#include "graph_simplification.hpp"
#include "coverage_handler.hpp"
#include "repeat_resolving.hpp"
#include "omni/omni_tools.hpp"
#include "seq_map.hpp"
#include "omni/ID_track_handler.hpp"
#include "omni/edges_position_handler.hpp"
#include "omni/total_labeler.hpp"
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "new_debruijn.hpp"
#include "config_struct.hpp"
#include "debruijn_stats.hpp"
#include "graphio.hpp"
#include "rectangleRepeatResolver.hpp"
#include "omni/distance_estimation.hpp"
#include "omni/loop_resolver.hpp"
#include "check_tools.hpp"
#include <cstdlib>

#include <boost/optional.hpp>
#include <boost/utility/in_place_factory.hpp>
#include <boost/utility/typed_in_place_factory.hpp>

//#include "dijkstra.hpp"

namespace debruijn_graph {

typedef io::IReader<io::SingleRead> SingleReadStream;
typedef io::IReader<io::PairedRead> PairedReadStream;
typedef io::ConvertingReaderWrapper UnitedStream;

using namespace omnigraph;

//For future refactoring
////---Graph Construction---
//void ConstructGraph();
//void SaveEnvAfterGC();
//void LoadEnvAfterGC();
//void CopyEnvAfterGC();
//
//void RunGraphConstruction(GraphPack gp) {
//	if (cfg::get().entry_point != debruijn::debruijn_config::after_construction) {
//		ConstructGraph(gp);
//		SaveEnvAfterGC(gp);
//	} else {
//		LoadEnvAfterGC(gp);
//		CopyEnvAfterGC(files);
//	}
//}
//
////---Paired Info Counting---
//void CountPairedInfo();
//void SaveEnvAfterPIC();
//void LoadEnvAfterPIC();
//void CopyEnvAfterPIC();
//
//void RunPairedInfoCounting() {
//	if (cfg::get().entry_point != debruijn::debruijn_config::after_pair_info_counting) {
//		RunGraphConstruction();
//		CountPairedInfo();
//		SaveEnvAfterPIC();
//	} else {
//		LoadEnvAfterPIC();
//		CopyEnvAfterPIC();
//	}
//}
//
//void Simplify();
//void SaveEnvAfterSimp();
//void LoadEnvAfterSimp();
//void CopyEnvAfterSimp();
//
//void RunSimplification() {
//	if (cfg::get().entry_point != debruijn::debruijn_config::after_simplification) {
//		RunPairedInfoCounting();
//		Simplify();
//		SaveEnvAfterSimp();
//	} else {
//		LoadEnvAfterSimp();
//		CopyEnvAfterSimp();
//	}
//}
//
//void RunSimplification();
//void EstimateDistances();
//void SaveEnvAfterDE();
//void LoadEnvAfterDE();
//void CopyEnvAfterDE();
//
//void RunDistanceEstimation() {
//	if (cfg::get().entry_point != debruijn::debruijn_config::after_distance_estimation) {
//		RunSimplification();
//		EstimateDistances();
//		SaveEnvAfterDE();
//	} else {
//		LoadEnvAfterDE();
//		CopyEnvAfterDE();
//	}
//}
//
//void ResolveRepeats();
//void SaveEnvAfterRR();
//void LoadEnvAfterRR();
//void CopyEnvAfterRR();
//
//void RunRepeatResolving() {
//	if (cfg::get().entry_point != debruijn::debruijn_config::after_repeat_resolving) {
//		RunDistanceEstimation();
//		ResolveRepeats();
//		SaveEnvAfterRR();
//	} else {
//		LoadEnvAfterRR();
//		CopyEnvAfterRR();
//	}
//}
//
//void RunAssembler() {
//	RunRepeatResolving();
//}

/*void FindZeros(PairedInfoIndex<Graph>& etalon_paired_index) {
 for (auto it = etalon_paired_index.begin(); it != etalon_paired_index.end(); ++it) {
 const vector<PairInfo<EdgeId>> infos = *it;
 for (auto it2 = infos.begin(); it2!=infos.end(); ++it2) {
 PairInfo<EdgeId> info = *it2;
 if (info.first == info.second && info.d == 0.) {
 cout << "FOUND ZEROS!!!" << endl;
 return;
 }
 }
 }
 }*/

template<size_t k, class Graph>
void SelectReadsForConsensus(Graph& etalon_graph, Graph& cur_graph,
		EdgeLabelHandler<Graph>& LabelsAfter,
		const EdgeIndex<k + 1, Graph>& index ,vector<SingleReadStream *>& reads
		, string& consensus_output_dir) {
	INFO("ReadMapping started");
	map<typename Graph::EdgeId, int> contigNumbers;
	int cur_num = 0;
	for (auto iter = cur_graph.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
		contigNumbers[*iter] = cur_num;
		cur_num++;
	}INFO(cur_num << "contigs");
	for (int i = 1; i < 3; i++) {
		int read_num = 0;
		osequencestream* mapped_reads[4000];
		for (int j = 0; j < cur_num; j++) {
			string output_filename = consensus_output_dir + ToString(j)
					+ "_reads" + ToString(i) + ".fa";
			osequencestream* tmp = new osequencestream(output_filename);
//			mapped_reads.push_back(tmp);
			mapped_reads[j] = tmp;
		}
		SingleReadMapper<k, Graph> rm(etalon_graph, index);
		while (!reads[i - 1]->eof()) {
			io::SingleRead cur_read;
			(*reads[i - 1]) >> cur_read;
			vector<typename Graph::EdgeId> res = rm.GetContainingEdges(
					cur_read);
			read_num++;
			TRACE(
					read_num<< " mapped to"<< res.size() <<" contigs :, read"<< cur_read.sequence());
//			map_quantity += res.size();
			for (size_t ii = 0; ii < res.size(); ii++) {
				TRACE("conting number "<< contigNumbers[res[ii]]);
				set<typename Graph::EdgeId> images =
						LabelsAfter.edge_inclusions[res[ii]];
				for (auto iter = images.begin(); iter != images.end(); ++iter)
					(*mapped_reads[contigNumbers[*iter]])
							<< cur_read.sequence();
			}
		}
	}
}

template<size_t k>
void CreateAndFillGraph(Graph& g, EdgeIndex<k + 1, Graph>& index
		, omnigraph::GraphLabeler<Graph>& labeler, PairedInfoIndex<Graph>& paired_index,
		PairedReadStream& stream, const Sequence& genome
		, EdgesPositionHandler<Graph> &EdgePos,
		PairedInfoIndex<Graph> &etalon_paired_index) {
	if (cfg::get().paired_mode) {
		if (cfg::get().etalon_info_mode) {
			ConstructGraphWithEtalonPairedInfo<k>(g, index,
					paired_index, stream, genome);
		} else {
			ConstructGraphWithPairedInfo<k>(g, index,
					paired_index, stream);
		}
		FillEtalonPairedIndex<k>(g, etalon_paired_index, index, genome);

	} else {
		UnitedStream united_stream(stream);
		ConstructGraphWithCoverage<k>(g, index,
				united_stream);
	}
//	ProduceInfo<k>(g, index, labeler, genome,
//			cfg::get().output_dir + "edge_graph.dot", "edge_graph");

//	FillEdgesPos<k>(g, index, genome, EdgePos);
}

template<size_t k>
void DeBruijnGraphTool(PairedReadStream& stream, const Sequence& genome,
		const string& work_tmp_dir, vector<io::Reader<io::SingleRead>*> reads) {

	using boost::optional;
	using boost::in_place;

	string output_folder = cfg::get().output_dir;

	INFO("Edge graph construction tool started");
	INFO("Start from " << cfg::get().start_from);
	INFO("Paired mode: " << (cfg::get().paired_mode ? "Yes" : "No"));
	INFO(
			"Etalon paired info mode: " << (cfg::get().etalon_info_mode ? "Yes" : "No"))
	mkdir(work_tmp_dir.c_str(),
			S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH | S_IWOTH);

	string graph_save_path = cfg::get().output_dir + "saves/";
	mkdir(graph_save_path.c_str(),
			S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH | S_IWOTH);

	Graph g(k);
	EdgeIndex<k + 1, Graph> index(g);
	IdTrackHandler<Graph> int_ids(g);
	EdgesPositionHandler<Graph> EdgePos(g);
	EdgesPosGraphLabeler<Graph> EdgePosLab(g, EdgePos);
	// if it's not paired_mode, then it'll be just unused variable -- takes O(1) to initialize from graph
	//	PairedInfoIndex<Graph> paired_index(g, 5);

	PairedInfoIndex<Graph> paired_index(g, 0);
	PairedInfoIndex<Graph> etalon_paired_index(g, 0);
	PairedInfoIndex<Graph> clustered_index(g);

	//experimental
//	KmerMapper<k + 1, Graph> kmer_mapper(g);
//	PairedInfoIndex<Graph> read_count_weight_paired_index(g, 0);
//	PairedInfoIndex<Graph> read_count_clustered_index(g);
	//experimental

//	int number_of_components = 0;
	bool graph_loaded = false;
	optional<TotalLabelerGraphStruct<Graph>> graph_struct;
	optional<TotalLabeler<Graph>> TotLab;

//	INFO(
//			"------Starting from: " << debruijn::debruijn_config::working_stage_name(cfg::get().entry_point) << "-----")

	if (cfg::get().start_from == "begin") {
		INFO("------Starting from Begin-----")

		graph_struct = in_place(boost::ref(g), &int_ids, &EdgePos);
		TotLab = in_place(&(*graph_struct));

		CreateAndFillGraph<k>(g, index, *TotLab, paired_index,
				stream, genome, EdgePos, etalon_paired_index);
		printGraph(g, int_ids, work_tmp_dir + "1_filled_graph", paired_index,
				EdgePos, &etalon_paired_index);
//		printGraph(g, int_ids, graph_save_path + "1_filled_graph", paired_index,
//				EdgePos, &etalon_paired_index);
		graph_loaded = true;


//		omnigraph::RealIdGraphLabeler<Graph> labeler(g, int_ids);

		omnigraph::WriteSimple(g, *TotLab, output_folder + "1_initial_graph.dot",
				"no_repeat_graph");
	}

	if (cfg::get().start_from == "after_filling") {
		scanConjugateGraph(&g, &int_ids, work_tmp_dir + "1_filled_graph",
				&paired_index, &EdgePos, &etalon_paired_index);
		graph_loaded = true;
		graph_struct = boost::in_place(boost::ref(g), &int_ids, &EdgePos);
		TotLab = in_place(&(*graph_struct));
	}

	if (cfg::get().start_from == "after_simplify") {
		scanConjugateGraph(&g, &int_ids, work_tmp_dir + "2_simplified_graph",
				&paired_index, &EdgePos, &etalon_paired_index,
				&clustered_index);
		graph_loaded = true;
		graph_struct = boost::in_place(boost::ref(g), &int_ids, &EdgePos);
		TotLab = in_place(&(*graph_struct));
	} else {
		if (graph_loaded) {
			SimplifyGraph<k>(g, index, *TotLab, 10, genome,
					output_folder/*, etalon_paired_index*/);
//			ProduceInfo<k>(g, index, *TotLab, genome,
//					output_folder + "simplified_graph.dot", "simplified_graph");

			//experimental
			//		FillPairedIndexWithReadCountMetric<k, ReadStream>(g, index, kmer_mapper, read_count_weight_paired_index, stream);
			//experimental

//			WriteGraphComponents<k>(g, index, *TotLab, genome,
//					output_folder + "graph_components" + "/", "graph.dot",
//					"graph_component", cfg::get().ds.IS);

//			number_of_components = PrintGraphComponents(
//					output_folder + "graph_components/graph", g,
//					cfg::get().ds.IS, int_ids, paired_index, EdgePos);
//
//
//			DistanceEstimator<Graph> estimator(g, paired_index,
//					cfg::get().ds.IS, cfg::get().ds.RL, cfg::get().de.delta,
//					cfg::get().de.linkage_distance, cfg::get().de.max_distance);
//			estimator.Estimate(clustered_index);
//
//			omnigraph::WriteSimple(g, *TotLab, output_folder + "2_simplified_graph.dot",
//					"no_repeat_graph");
//
//			//todo think if we need this save
//			printGraph(g, int_ids, graph_save_path + "repeats_resolved_before",
//					paired_index, EdgePos/*, &read_count_weight_paired_index*/);
//
			printGraph(g, int_ids, work_tmp_dir + "2_simplified_graph",
					paired_index, EdgePos, &etalon_paired_index,
					&clustered_index/*, &read_count_weight_paired_index*/);
//			printGraph(g, int_ids, output_folder + "2_simplified_graph",
//					clustered_index, EdgePos, &etalon_paired_index,
//					&clustered_index/*, &read_count_weight_paired_index*/);
		}
	}

//	if (cfg::get().paired_mode) {
//		if (graph_loaded) {
//			CountPairedInfoStats(g, paired_index, etalon_paired_index,
//					output_folder);
//
//			//experiment
//			//			DistanceEstimator<Graph> estimator2(g, read_count_weight_paired_index, insert_size,
//			//					max_read_length, cfg::get().de.delta,
//			//					cfg::get().de.linkage_distance,
//			//					cfg::get().de.max_distance);
//			//			estimator2.Estimate(read_count_clustered_index);
//			//experiment
//
//			CountClusteredPairedInfoStats(g, paired_index, clustered_index,
//					etalon_paired_index, output_folder);
//
//		}
//
//// after_simplify
//		INFO("before ResolveRepeats");
//
//		NCGraph new_graph(k);
//		IdTrackHandler<NCGraph> NewIntIds(new_graph, int_ids.MaxVertexId(),
//				int_ids.MaxEdgeId());
//		PairedInfoIndex<NCGraph> new_index(new_graph);
//		EdgeIndex<k + 1, NCGraph> new_edge_index(new_graph);
//		EdgesPositionHandler<NCGraph> EdgePosBefore(new_graph);
//
//		scanNCGraph(new_graph, NewIntIds, work_tmp_dir + "2_simplified_graph",
//				(PairedInfoIndex<NCGraph>*) 0, EdgePosBefore,
//				(PairedInfoIndex<NCGraph>*) 0, &new_index);
//
//		if (cfg::get().start_from == "after_simplify"
//				|| cfg::get().start_from == "before_resolve") {
//			//todo ask Shurik if graph is not empty here
////			WriteGraphComponents<k>(g, index, genome,
////					output_folder + "graph_components" + "/", "graph.dot",
////					"graph_component", cfg::get().ds.IS);
//
//			mkdir((output_folder + "graph_components" + "/").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH | S_IWOTH);
//		}
//
//
//		number_of_components = PrintGraphComponents(
//				output_folder + "graph_components/graphCl", new_graph,
//				cfg::get().ds.IS, NewIntIds, new_index, EdgePosBefore);
//
//		RealIdGraphLabeler<NCGraph> IdTrackLabelerAfter(new_graph, NewIntIds);
//
//		INFO("repeat resolved graph written");
//
//		NonconjugateDeBruijnGraph resolved_graph(k);
//		IdTrackHandler<NCGraph> Resolved_IntIds(resolved_graph);
//		EdgesPositionHandler<NCGraph> EdgePosAfter(resolved_graph);
//		EdgeLabelHandler<NCGraph> LabelsAfter(resolved_graph, new_graph);
//
//		DEBUG("New index size: "<< new_index.size());
//		if (cfg::get().rectangle_mode) {
//			void RectangleResolve(
//					PairedInfoIndex<NonconjugateDeBruijnGraph>& index,
//					NonconjugateDeBruijnGraph& graph,
//					const string& work_tmp_dir, const string& output_folder);
//			RectangleResolve(new_index, new_graph, work_tmp_dir, output_folder);
//		}
//
//		ResolveRepeats(new_graph, NewIntIds, new_index, EdgePosBefore,
//				resolved_graph, Resolved_IntIds, EdgePosAfter,
//				output_folder + "resolve/", LabelsAfter);
//
//		INFO("Total labeler start");
//		TotalLabelerGraphStruct<NCGraph> graph_struct_before(new_graph,
//				&NewIntIds, &EdgePosBefore, NULL);
//		TotalLabelerGraphStruct<NCGraph> graph_struct_after(resolved_graph,
//				&Resolved_IntIds, &EdgePosAfter, &LabelsAfter);
//		TotalLabeler<NCGraph> TotLabAfter(&graph_struct_after,
//				&graph_struct_before);
//
//		omnigraph::WriteSimple(resolved_graph, TotLabAfter, output_folder + "3_resolved_graph.dot",
//				"no_repeat_graph");
//
//		INFO("Total labeler finished");
//
//		RealIdGraphLabeler<NCGraph> IdTrackLabelerResolved(resolved_graph,
//				Resolved_IntIds);
//
//		EdgesPosGraphLabeler<NCGraph> EdgePosLAfterLab(resolved_graph,
//				EdgePosAfter);
//
//		EdgesLabelsGraphLabeler<NCGraph> LabelLabler(resolved_graph,
//				LabelsAfter);
//		INFO("---Clearing resolved graph---");
//		for (int i = 0; i < 3; i++) {
//			ClipTipsForResolver(resolved_graph);
//			RemoveBulges2(resolved_graph);
//			RemoveLowCoverageEdgesForResolver(resolved_graph);
//		}INFO("---Cleared---");
//		INFO("---Output Contigs---");
//		OutputContigs(resolved_graph,
//				output_folder + "contigs_before_enlarge.fasta");
//
//		omnigraph::WriteSimple(resolved_graph, TotLabAfter, output_folder + "4_cleared_graph.dot",
//				"no_repeat_graph");
//
//		one_many_contigs_enlarger<NCGraph> N50enlarger(resolved_graph);
//		N50enlarger.Loops_resolve();
//
//		omnigraph::WriteSimple(resolved_graph, TotLabAfter, output_folder + "5_unlooped_graph.dot",
//				"no_repeat_graph");
//		OutputContigs(resolved_graph, output_folder + "contigs_unlooped.fasta");
//
//		N50enlarger.one_many_resolve_with_vertex_split();
//
//		omnigraph::WriteSimple(resolved_graph, TotLabAfter, output_folder + "6_finished_graph.dot",
//				"no_repeat_graph");
//
//		OutputContigs(resolved_graph, output_folder + "contigs_final.fasta");
//		string consensus_folder = output_folder + "consensus/";
//
////		OutputSingleFileContigs(resolved_graph, consensus_folder);
////		SelectReadsForConsensus<k, NCGraph>(new_graph, resolved_graph, LabelsAfter, new_edge_index, reads, consensus_folder);
//
//		OutputContigs(new_graph,
//				output_folder + "contigs_before_resolve.fasta");
//
//		if (number_of_components > 0) {
//			string output_comp = output_folder + "resolved_comp";
//			mkdir(output_comp.c_str(),
//					S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH | S_IWOTH);
//
//			for (int i = 1; i <= number_of_components; i++)
//				ResolveOneComponent(output_folder + "graph_components/",
//						output_comp + "/", i, k);
//		}
//	}
//	if (!cfg::get().paired_mode)
//		OutputContigs(g, output_folder + "contigs.fasta");
	INFO("Tool finished");
}
}

#endif /* LAUNCH_HPP_ */
