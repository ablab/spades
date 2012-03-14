#pragma once

#include "repeat_resolver.hpp"
#include "graphio.hpp"
#include "omni/one_many_contigs_enlarger.hpp"
#include "omni/loop_resolver.hpp"
#include <sys/types.h>
#include <sys/stat.h>

namespace debruijn_graph {

template<class Graph>
void ResolveRepeats(Graph &g, IdTrackHandler<Graph> &old_IDs,
		PairedInfoIndex<Graph> &info, EdgesPositionHandler<Graph> &edges_pos,
		Graph &new_graph, IdTrackHandler<Graph> &new_IDs,
		EdgesPositionHandler<Graph> &edges_pos_new, const string& output_folder,
		EdgeLabelHandler<Graph> &LabelsAfter) {
	INFO("SUBSTAGE == Resolving primitive repeats");
	for (auto e_iter = g.SmartEdgeBegin(); !e_iter.IsEnd();
			++e_iter) {
		{
			if ( (g.EdgeStart(*e_iter) == g.EdgeEnd(*e_iter))){
				if ((g.length(*e_iter) > 1) ){
					size_t half_len = g.length(*e_iter) / 2;
					g.SplitEdge(*e_iter, half_len);
				}
				else {
					WARN("Loop of length 1 detected");
				}
			}
		}
	}
	DeletedVertexHandler<Graph> tmp_deleted_handler(new_graph);
	TRACE("deleted handler created");
	RepeatResolver<Graph> repeat_resolver(g, old_IDs, info, edges_pos,
			new_graph, new_IDs, edges_pos_new, tmp_deleted_handler, LabelsAfter);
	make_dir(output_folder);
	auto edge_labels = repeat_resolver.GetEdgeLabels();
	LabelsAfter.FillLabels(edge_labels);

	repeat_resolver.ResolveRepeats(output_folder);
	DEBUG("Primitive repeats resolved");
}

//void ResolveOneComponent(const string& load_from_dir, const string& save_to_dir,
//		int component_id, int k) {
//	string load_from = ConstructComponentName(load_from_dir + "/graphCl",
//			component_id);
//	string save_to = ConstructComponentName(save_to_dir + "/graph",
//			component_id);
//
//	string save_resolving_history = ConstructComponentName(
//			save_to_dir + "/resolve", component_id);
//	make_dir(save_resolving_history);
//
//	NCGraph new_graph(k);
//	IdTrackHandler<NCGraph> NewIntIds(new_graph);
//	PairedInfoIndex<NCGraph> new_index(new_graph);
//	EdgesPositionHandler<NCGraph> EdgePosBefore(new_graph);
//	scanNCGraph(new_graph, NewIntIds, load_from, &new_index, EdgePosBefore);
//
//	RealIdGraphLabeler<NCGraph> IdTrackLabelerAfter(new_graph, NewIntIds);
//
//	omnigraph::WriteSimple(new_graph, IdTrackLabelerAfter, save_to + "_before.dot", "no_repeat_graph");
//
//	NonconjugateDeBruijnGraph resolved_graph(k);
//	IdTrackHandler<NCGraph> Resolved_IntIds(resolved_graph);
//	EdgesPositionHandler<NCGraph> EdgePosAfter(resolved_graph);
//	EdgeLabelHandler<NCGraph> LabelsAfter(resolved_graph, new_graph);
//
//	ResolveRepeats(new_graph, NewIntIds, new_index, EdgePosBefore,
//			resolved_graph, Resolved_IntIds, EdgePosAfter,
//			save_resolving_history + "/", LabelsAfter);
//
//	RealIdGraphLabeler<NCGraph> IdTrackLabelerResolved(resolved_graph,
//			Resolved_IntIds);
//
//	TotalLabelerGraphStruct<NCGraph> graph_struct_before(new_graph,
//				&NewIntIds, &EdgePosBefore, NULL);
//	TotalLabelerGraphStruct<NCGraph> graph_struct_after(resolved_graph,
//				&Resolved_IntIds, &EdgePosAfter, &LabelsAfter);
//	TotalLabeler<NCGraph> TotLabAfter(&graph_struct_after,
//				&graph_struct_before);
//
//	omnigraph::WriteSimple(resolved_graph, TotLabAfter, save_to + "_after.dot", "no_repeat_graph");
//
//	EdgesPosGraphLabeler<NCGraph> EdgePosLAfterLab(resolved_graph,
//			EdgePosAfter);
//
//	omnigraph::WriteSimple(resolved_graph, TotLabAfter,
//			save_resolving_history + "/repeats_resolved_after_pos.dot",
//			"no_repeat_graph");
//
//	ClipTips(resolved_graph);
//	RemoveLowCoverageEdgesForResolver(resolved_graph);
//
//	omnigraph::WriteSimple(resolved_graph, TotLabAfter,
//			save_resolving_history
//					+ "/repeats_resolved_after_und_cleared_pos.dot",
//			"no_repeat_graph");
//	omnigraph::WriteSimple(resolved_graph, TotLabAfter,
//			save_resolving_history + "/repeats_resolved_und_cleared.dot",
//			"no_repeat_graph");
//	one_many_contigs_enlarger<NCGraph> N50enlarger(resolved_graph, cfg::get().ds.IS);
//	N50enlarger.one_many_resolve_with_vertex_split();
//	omnigraph::WriteSimple(resolved_graph, IdTrackLabelerResolved, save_to + "_finished.dot", "no_repeat_graph");
//}

//void ConjugateResolveOneComponent(const string& load_from_dir, const string& save_to_dir,
//		int component_id, int k) {
//	string file_name = ConstructComponentName(load_from_dir + "/graphCl",
//			component_id);
//	string save_to = ConstructComponentName(save_to_dir + "/graph",
//			component_id);
//
//	string save_resolving_history = ConstructComponentName(
//			save_to_dir + "/resolve", component_id);
//	make_dir(save_resolving_history);
//
//	Graph graph(k);
//	IdTrackHandler<Graph> int_ids(graph);
//	PairedInfoIndex<Graph> paired_index(graph);
//	EdgesPositionHandler<Graph> edge_pos(graph, cfg::get().pos.max_single_gap);
//
//	ConjugateDataScanner<Graph> scanner(graph, int_ids);
//	ScanBasicGraph(file_name, scanner);
//	scanner.loadPositions(file_name, edge_pos);
//	scanner.loadPaired(file_name, paired_index);
//
//	RealIdGraphLabeler<Graph> id_track_labeler(graph, int_ids);
//
//	omnigraph::WriteSimple(graph, id_track_labeler, save_to + "_before.dot", "no_repeat_graph");
//
//	ConjugateDeBruijnGraph resolved_graph(k);
//	IdTrackHandler<Graph> resolved_int_ids(resolved_graph);
//	EdgesPositionHandler<Graph> resolved_edge_pos(resolved_graph, cfg::get().pos.max_single_gap);
//	EdgeLabelHandler<Graph> resolve_mapper(resolved_graph, graph);
//
//	ResolveRepeats(graph, int_ids, paired_index, edge_pos,
//			resolved_graph, resolved_int_ids, resolved_edge_pos,
//			save_resolving_history + "/", resolve_mapper);
//
//	RealIdGraphLabeler<Graph> resolved_id_track_labeler(resolved_graph,
//			resolved_int_ids);
//	omnigraph::WriteSimple(resolved_graph, resolved_id_track_labeler, save_to + "_after.dot", "no_repeat_graph");
//
//	EdgesPosGraphLabeler<Graph> resolved_edge_pos_labeler(resolved_graph,
//			resolved_edge_pos);
//
//	omnigraph::WriteSimple(resolved_graph, resolved_edge_pos_labeler,
//			save_resolving_history + "/repeats_resolved_after_pos.dot",
//			"no_repeat_graph");
//
//	ClipTips(resolved_graph);
//	RemoveLowCoverageEdgesForResolver(resolved_graph);
//
//	omnigraph::WriteSimple(resolved_graph, resolved_edge_pos_labeler,
//			save_resolving_history
//					+ "/repeats_resolved_after_und_cleared_pos.dot",
//			"no_repeat_graph");
////	omnigraph::WriteSimple(resolved_graph, resolved_id_track_labeler,
////			save_resolving_history + "/repeats_resolved_und_cleared.dot",
////			"no_repeat_graph");
////	one_many_contigs_enlarger<Graph> N50enlarger(resolved_graph, cfg::get().ds.IS);
////	N50enlarger.one_many_resolve_with_vertex_split();
////	omnigraph::WriteSimple(resolved_graph, resolved_id_track_labeler, save_to + "_finished.dot", "no_repeat_graph");
//}


//void RectangleResolve(PairedInfoIndex<NonconjugateDeBruijnGraph>& index,
//		NonconjugateDeBruijnGraph& graph, const string& work_tmp_dir,
//		const string& output_folder) {
//
//	NonconjugateDeBruijnGraph resolvedGraph(graph.k());
//	typedef NonconjugateDeBruijnGraph::EdgeId NCEdgeId;
//	PairInfoIndexData<NCEdgeId> piid;
//	for (auto iter = index.begin(); iter != index.end(); ++iter) {
//		vector<PairInfo<NCEdgeId> > pi = *iter;
//		for (size_t i = 0; i < pi.size(); ++i) {
//			if (pi[i].d >= 0)
//				piid.AddPairInfo(pi[i], 1);
//		}
//	}
//	RectangleRepeatResolver<NonconjugateDeBruijnGraph> rectangleResolver(graph,
//			piid, resolvedGraph, (size_t) 30);
//	rectangleResolver.Process();
//
//	IdTrackHandler<NCGraph> Resolved_IntIds(resolvedGraph);
//	RealIdGraphLabeler<NCGraph> IdTrackLabelerResolved(resolvedGraph,
//			Resolved_IntIds);
//
//	omnigraph::WriteSimple(resolvedGraph, IdTrackLabelerResolved, work_tmp_dir + "rectgraph.dot", "rectgraph");
//	INFO("rect graph written: " + work_tmp_dir + "rectgraph.dot");
//
//	for (auto iter = resolvedGraph.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
//		INFO("COV:" << resolvedGraph.coverage(*iter));
//	}
//
////	ClipTips(resolvedGraph);
////	RemoveLowCoverageEdges(resolvedGraph);
////
////	ClipTips(resolvedGraph);
////	RemoveLowCoverageEdges(resolvedGraph);
////	see if two methods result in the same graph.
//
//	for (int i = 0; i < 3; i++) {
//		ClipTips(resolvedGraph);
//		RemoveBulges2(resolvedGraph);
//		RemoveLowCoverageEdgesForResolver(resolvedGraph);
//
//	}
//	LoopResolver<NCGraph> loopResolver(resolvedGraph, 0.5);
//	loopResolver.ResolveLoops();
//
//	one_many_contigs_enlarger<NCGraph> N50enlarger(resolvedGraph, cfg::get().ds.IS);
//	N50enlarger.one_many_resolve_with_vertex_split();
//	N50enlarger.Loops_resolve();
//	omnigraph::Compressor<NCGraph> compressor(resolvedGraph);
//	compressor.CompressAllVertices();
//	omnigraph::Cleaner<NCGraph> cleaner(resolvedGraph);
//	cleaner.Clean();
//
//	IdTrackHandler<NCGraph> idTrackerAfter(resolvedGraph);
//	RealIdGraphLabeler<NCGraph> idLabelAfter(resolvedGraph, idTrackerAfter);
//
//	omnigraph::WriteSimple(resolvedGraph, idLabelAfter, work_tmp_dir + "rectgraphAfter.dot",
//			"rectgraphAfter");
//	INFO("rect graph written: " + work_tmp_dir + "rectgraphAfter.dot");
//
//	EmptyGraphLabeler<NonconjugateDeBruijnGraph> emptyLabeler;
//	omnigraph::WriteSimple(graph, emptyLabeler, work_tmp_dir + "beforerectgraph.dot",
//			"beforerectgraph");
//	INFO("rect graph written: " + work_tmp_dir + "beforerectgraph.dot");
//
//	OutputContigs(resolvedGraph, output_folder + "rectcontig.fasta");
//	OutputContigs(graph, output_folder + "before-rectcontig.fasta");
//}

}
