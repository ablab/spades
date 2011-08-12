/*
 * launch.hpp
 *
 *  Created on: May 6, 2011
 *      Author: sergey
 */

#ifndef LAUNCH_HPP_
#define LAUNCH_HPP_

#include "io/reader.hpp"
#include "io/rc_reader_wrapper.hpp"
#include "io/cutting_reader_wrapper.hpp"
#include "io/converting_reader_wrapper.hpp"
#include "visualization_utils.hpp"

//#include "debruijn_graph.hpp"
#include "edge_labels_handler.hpp"
#include "paired_info.hpp"
#include "graph_construction.hpp"
#include "graph_simplification.hpp"
#include "coverage_handler.hpp"
#include "repeat_resolving.hpp"
#include "omni_tools.hpp"
#include "seq_map.hpp"
#include "ID_track_handler.hpp"
#include "edges_position_handler.hpp"
#include "total_labeler.hpp"
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "new_debruijn.hpp"
#include "config.hpp"
#include "debruijn_stats.hpp"
#include "graphio.hpp"
#include "rectangleRepeatResolver.hpp"
#include "distance_estimation.hpp"
#include "loop_resolver.hpp"
#include "check_tools.hpp"
#include <cstdlib>
//#include "dijkstra.hpp"

namespace debruijn_graph {

typedef io::Reader<io::SingleRead> SingleReadStream;
using namespace omnigraph;


template<size_t k, class Graph>
void SelectReadsForConsensus(Graph& g, const EdgeIndex<k + 1, Graph>& index ,vector<SingleReadStream *>& reads, string& consensus_output_dir){
	INFO("ReadMapping started");
	map<typename Graph::EdgeId, int> contigNumbers;
	int cur_num = 0;
	for (auto iter = g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
		contigNumbers[*iter] = cur_num;
		cur_num++;
	}
	INFO(cur_num << "contigs");
	for( int i = 1; i < 3; i ++ ) {
		int read_num = 0;
		osequencestream* mapped_reads[2000];
		for(int j = 0; j < cur_num; j++){
			string output_filename = consensus_output_dir + ToString(j) + "_reads" + ToString(i)+".fa";
			osequencestream* tmp = new osequencestream(output_filename);
//			mapped_reads.push_back(tmp);
			mapped_reads[j] = tmp;
		}
		SingleReadMapper<k, Graph> rm(g, index);
		while (!reads[i - 1]->eof()) {
			io::SingleRead cur_read;
			(*reads[i - 1]) >> cur_read;
			vector<typename Graph::EdgeId> res = rm.GetContainingEdges(cur_read);
			read_num++;
			TRACE(read_num<< " mapped to"<< res.size() <<" contigs :, read"<< cur_read.sequence());
//			map_quantity += res.size();
			for(size_t ii = 0; ii < res.size(); ii++) {
				TRACE("conting number "<< contigNumbers[res[ii]]);
				(*mapped_reads[contigNumbers[res[ii]]]) << cur_read.sequence();
			}
		}
	}
}

template<size_t k, class ReadStream>
void DeBruijnGraphTool(ReadStream& stream, const Sequence& genome,
		bool paired_mode, bool rectangle_mode, bool etalon_info_mode,
		bool from_saved, size_t insert_size, size_t max_read_length,
		const string& output_folder, const string& work_tmp_dir, vector<SingleReadStream* > reads) {

	INFO("Edge graph construction tool started");
	INFO("Paired mode: " << (paired_mode ? "Yes" : "No"));
	INFO("Etalon paired info mode: " << (etalon_info_mode ? "Yes" : "No"))INFO(
			"From file: " << (from_saved ? "Yes" : "No"))
	mkdir(work_tmp_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH | S_IWOTH);

	string graph_save_path = output_folder+"saves/";
	mkdir(graph_save_path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH | S_IWOTH);

	Graph g(k);
	EdgeIndex<k + 1, Graph> index(g);
	IdTrackHandler<Graph> IntIds(g);
	EdgesPositionHandler<Graph> EdgePos(g);
	EdgesPosGraphLabeler<Graph> EdgePosLab(g, EdgePos);
	// if it's not paired_mode, then it'll be just unused variable -- takes O(1) to initialize from graph
	//	PairedInfoIndex<Graph> paired_index(g, 5);
	PairedInfoIndex<Graph> paired_index(g, 0);
	PairedInfoIndex<Graph> etalon_paired_index(g, 0);

	PairedInfoIndex<Graph> clustered_index(g);
	int number_of_components = 0;

	if (!from_saved) {
		if (paired_mode) {
			if (etalon_info_mode) {
				ConstructGraphWithEtalonPairedInfo<k, ReadStream> (g, index,
						paired_index, stream, insert_size, max_read_length,
						genome);
			} else {
				ConstructGraphWithPairedInfo<k, ReadStream> (g, index,
						paired_index, stream);
			}
			FillEtalonPairedIndex<k> (g, etalon_paired_index, index,
					insert_size, max_read_length, genome);
		} else {
			typedef io::ConvertingReaderWrapper UnitedStream;
			UnitedStream united_stream(&stream);
			ConstructGraphWithCoverage<k, UnitedStream> (g, index,
					united_stream);
		}

		ProduceInfo<k> (g, index, genome, output_folder + "edge_graph.dot",
				"edge_graph");
		FillEdgesPos<k> (g, index, genome, EdgePos);

		TotalLabelerGraphStruct<Graph> graph_struct(g, &IntIds, &EdgePos, NULL);
		TotalLabeler<Graph> TotLab(&graph_struct, NULL);


		omnigraph::WriteSimple(output_folder + "1_initial_graph.dot",
				"no_repeat_graph", g, TotLab);

//		omnigraph::WriteSimple(output_folder + "before_simplification_pos.dot",
//				"no_repeat_graph", g, EdgePosLab);

		printGraph(g, IntIds, graph_save_path + "first_graph", paired_index,
				EdgePos);

		SimplifyGraph<k> (g, index, 3, genome, output_folder/*, etalon_paired_index*/);
		ProduceInfo<k> (g, index, genome,
				output_folder + "simplified_graph.dot", "simplified_graph");

		WriteGraphComponents<k> (g, index, genome,
				output_folder + "graph_components" + "/", "graph.dot",
				"graph_component", insert_size);

		number_of_components = PrintGraphComponents(
				output_folder + "graph_components/graph", g, insert_size,
				IntIds, paired_index, EdgePos);

		if (paired_mode) {
			CountPairedInfoStats(g, insert_size, max_read_length, paired_index,
					etalon_paired_index, output_folder);
		}

//		omnigraph::WriteSimple(
//				output_folder + "repeats_resolved_before_poslab.dot",
//				"no_repeat_graph", g, EdgePosLab);
//		omnigraph::WriteSimple(
//				work_tmp_dir + "repeats_resolved_before_poslab.dot",
//				"no_repeat_graph", g, EdgePosLab);

		printGraph(g, IntIds, graph_save_path + "repeats_resolved_before",
				paired_index, EdgePos);

		omnigraph::WriteSimple(output_folder + "2_simplified_graph.dot",
				"no_repeat_graph", g, TotLab);

		if (paired_mode) {
			DistanceEstimator<Graph> estimator(g, paired_index, insert_size,
					max_read_length, CONFIG.read<size_t> ("de_delta"),
					CONFIG.read<size_t> ("de_linkage_distance"),
					CONFIG.read<size_t> ("de_max_distance"));
			estimator.Estimate(clustered_index);
			CountClusteredPairedInfoStats(g, insert_size, max_read_length,
					paired_index, clustered_index, etalon_paired_index,
					output_folder);
		}

		number_of_components = PrintGraphComponents(
				output_folder + "graph_components/graphCl", g, insert_size,
				IntIds, clustered_index, EdgePos);

	}
	//	if (paired_mode) {
	//		paired_index.OutputData(output_folder + "edges_dist.txt");
	//
	//		SimpleOfflineClusterer<Graph> clusterer(paired_index);
	//		PairedInfoIndex<Graph> clustered_paired_index(g);
	//		clusterer.cluster(clustered_paired_index);
	//	}


	if (paired_mode) {
		if (!from_saved) {
			INFO("before ResolveRepeats");
			RealIdGraphLabeler<Graph> IdTrackLabelerBefore(g, IntIds);
//			omnigraph::WriteSimple(
//					output_folder + "repeats_resolved_before.dot",
//					"no_repeat_graph", g, IdTrackLabelerBefore);
			printGraph(g, IntIds, work_tmp_dir + "graph", clustered_index,
					EdgePos);
			printGraph(g, IntIds, output_folder + "graph", clustered_index,
					EdgePos);
		}

		NCGraph new_graph(k);
		IdTrackHandler<NCGraph> NewIntIds(new_graph, IntIds.MaxVertexId(),
				IntIds.MaxEdgeId());
		PairedInfoIndex<NCGraph> new_index(new_graph);
		EdgeIndex<k+1, NCGraph> new_edge_index(new_graph);
		EdgesPositionHandler<NCGraph> EdgePosBefore(new_graph);

		/*Graph conj_copy_graph(k);
		IdTrackHandler<Graph> conj_IntIds(conj_copy_graph,
				IntIds.MaxVertexId(), IntIds.MaxEdgeId());
		PairedInfoIndex<Graph> conj_copy_index(conj_copy_graph);

		 scanConjugateGraph(conj_copy_graph, conj_IntIds,
		 work_tmp_dir + "graph", conj_copy_index);
		 printGraph(conj_copy_graph, conj_IntIds, work_tmp_dir + "graph_copy",
		 conj_copy_index);
		 */
		scanNCGraph(new_graph, NewIntIds, work_tmp_dir + "graph", new_index,
				EdgePosBefore);

		RealIdGraphLabeler<NCGraph> IdTrackLabelerAfter(new_graph, NewIntIds);

//		omnigraph::WriteSimple(work_tmp_dir + "repeats_resolved_nonconjugate_copy.dot",
//				"no_repeat_graph", new_graph, IdTrackLabelerAfter);
		INFO("repeat resolved graph written");

		NonconjugateDeBruijnGraph resolved_graph(k);
		IdTrackHandler<NCGraph> Resolved_IntIds(resolved_graph);
		EdgesPositionHandler<NCGraph> EdgePosAfter(resolved_graph);
		EdgeLabelHandler<NCGraph> LabelsAfter(resolved_graph, new_graph);

		DEBUG("New index size: "<< new_index.size());
		if (rectangle_mode) {
			void RectangleResolve(
					PairedInfoIndex<NonconjugateDeBruijnGraph>& index,
					NonconjugateDeBruijnGraph& graph,
					const string& work_tmp_dir, const string& output_folder);
			RectangleResolve(new_index, new_graph, work_tmp_dir, output_folder);
		}

		ResolveRepeats(new_graph, NewIntIds, new_index, EdgePosBefore,
				resolved_graph, Resolved_IntIds, EdgePosAfter,
				output_folder + "resolve/", LabelsAfter);

		INFO("Total labeler start");
		TotalLabelerGraphStruct<NCGraph> graph_struct_before(new_graph, &NewIntIds, &EdgePosBefore, NULL);
		TotalLabelerGraphStruct<NCGraph> graph_struct_after(resolved_graph, &Resolved_IntIds, &EdgePosAfter, &LabelsAfter);
		TotalLabeler<NCGraph> TotLabAfter(&graph_struct_after, &graph_struct_before);

		omnigraph::WriteSimple(output_folder + "3_resolved_graph.dot",
				"no_repeat_graph", resolved_graph, TotLabAfter);

		INFO("Total labeler finished");

		RealIdGraphLabeler<NCGraph> IdTrackLabelerResolved(resolved_graph,
				Resolved_IntIds);

//		omnigraph::WriteSimple(work_tmp_dir + "repeats_resolved_after.dot",
//				"no_repeat_graph", resolved_graph, IdTrackLabelerResolved);
//		omnigraph::WriteSimple(output_folder + "repeats_resolved_after.dot",
//				"no_repeat_graph", resolved_graph, IdTrackLabelerResolved);

		EdgesPosGraphLabeler<NCGraph> EdgePosLAfterLab(resolved_graph,
				EdgePosAfter);

//		omnigraph::WriteSimple(work_tmp_dir + "repeats_resolved_after_pos.dot",
//		"no_repeat_graph", resolved_graph, EdgePosLAfterLab);
//		omnigraph::WriteSimple(
//				output_folder + "repeats_resolved_after_pos.dot",
//				"no_repeat_graph", resolved_graph, EdgePosLAfterLab);

		EdgesLabelsGraphLabeler<NCGraph> LabelLabler(resolved_graph, LabelsAfter);

//		omnigraph::WriteSimple(
//				output_folder + "resolved_labels_1.dot",
//				"no_repeat_graph", resolved_graph, LabelLabler);

		for(int i = 0; i < 3; i ++) {
			ClipTipsForResolver(resolved_graph);
			RemoveBulges2(resolved_graph);
			RemoveLowCoverageEdgesForResolver(resolved_graph);
		}
//		omnigraph::WriteSimple(
//				output_folder + "resolved_labels_2.dot",
//				"no_repeat_graph", resolved_graph, LabelLabler);
//		omnigraph::WriteSimple(
//				work_tmp_dir + "repeats_resolved_after_und_cleared_pos.dot",
//				"no_repeat_graph", resolved_graph, EdgePosLAfterLab);
//		omnigraph::WriteSimple(
//				output_folder + "repeats_resolved_after_und_cleared_pos.dot",
//				"no_repeat_graph", resolved_graph, EdgePosLAfterLab);

//		omnigraph::WriteSimple(
//				work_tmp_dir + "repeats_resolved_und_cleared.dot",
//				"no_repeat_graph", resolved_graph, IdTrackLabelerResolved);
//		omnigraph::WriteSimple(
//				output_folder + "repeats_resolved_und_cleared.dot",
//				"no_repeat_graph", resolved_graph, IdTrackLabelerResolved);

		OutputContigs(resolved_graph, output_folder + "contigs_before_enlarge.fasta");

		omnigraph::WriteSimple(output_folder + "4_cleared_graph.dot",
				"no_repeat_graph", resolved_graph, TotLabAfter);

//		omnigraph::WriteSimple(
//				output_folder + "repeats_resolved_und_cleared.dot",
//				"no_repeat_graph", resolved_graph, IdTrackLabelerResolved);

		one_many_contigs_enlarger<NCGraph> N50enlarger(resolved_graph);
		N50enlarger.one_many_resolve_with_vertex_split();

		omnigraph::WriteSimple(output_folder + "4_finished_graph.dot",
				"no_repeat_graph", resolved_graph, TotLabAfter);

//		omnigraph::WriteSimple(
//				output_folder + "resolved_labels_3.dot",
//				"no_repeat_graph", resolved_graph, LabelLabler);
//
//
//		omnigraph::WriteSimple(
//				output_folder
//						+ "repeats_resolved_und_cleared_und_simplified.dot",
//				"no_repeat_graph", resolved_graph, IdTrackLabelerResolved);
//		INFO("repeat resolved grpah written");

		//		SimplifyGRaph<k>(resolved_graph, aux_index, 3, genome, output_folder);

		//CountStats<k, NCGraph> (resolved_graph, aux_index, genome);
		//		ProduceNonconjugateInfo<k> (resolved_graph, aux_index, genome, output_folder + "repeats_resolved.dot",
		//				"no_repeat_graph");
		//		ProduceNonconjugateInfo<k> (resolved_graph, aux_index, genome, work_tmp_dir + "repeats_resolved.dot",
		//				"no_repeat_graph");sss

		OutputContigs(resolved_graph, output_folder + "contigs_final.fasta");
		string consensus_folder = output_folder + "consensus/";

//		OutputSingleFileContigs(new_graph, consensus_folder);
//		SelectReadsForConsensus<k, NCGraph>(new_graph, new_edge_index, reads, consensus_folder);

		OutputContigs(new_graph, output_folder + "contigs_before_resolve.fasta");

		if (number_of_components > 0) {
			string output_comp = output_folder + "resolved_comp";
			mkdir(output_comp.c_str(),
					S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH | S_IWOTH);

			for (int i = 1; i <= number_of_components; i++)
				ResolveOneComponent(output_folder + "graph_components/",
						output_comp + "/", i, k);
		}

	}
	if (!paired_mode)
		OutputContigs(g, output_folder + "contigs.fasta");
	INFO("Tool finished");
}

}

#endif /* LAUNCH_HPP_ */
