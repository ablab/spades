#pragma once
#include "omni/matepair_transformer.hpp"
#include "omni/path_set.hpp"
#include "new_debruijn.hpp"
#include "graph_pack.hpp"
#include "utils.hpp"
#include "omni/omni_utils.hpp"

#include "omni/omni_tools.hpp"
#include "omni/omnigraph.hpp"

#include "omni/edges_position_handler.hpp"
#include "omni/total_labeler.hpp"
#include "path_set_stats.hpp"

namespace debruijn_graph{
template <class Graph>

class PathSetGraphConstructor {

typedef typename Graph::EdgeId EdgeId;
typedef typename Graph::VertexId VertexId;

typedef vector<EdgeId > Path;
const Graph& g_;
const PairedInfoIndex<Graph>& pair_info_;

public:
PathSetGraphConstructor(const Graph& g,const PairedInfoIndex<Graph>& pair_info, Graph& new_graph, IdTrackHandler<Graph>& newIds, TotalLabeler<Graph>& tot_labeler_after): g_(g), pair_info_(pair_info) {
	PathSetIndexData<EdgeId> PII ;
	PathSetIndexData<EdgeId> PIIFilter ;
	MatePairTransformer<Graph> transformer(g_, pair_info_);
	transformer.Transform(PII);
	PathSetIndex<EdgeId> PI(PII);
	PI.RemovePrefixes(PIIFilter);


    PathSetStats<Graph>  pathsetStatistic(g,PII , PIIFilter);
    pathsetStatistic.Count();


	for(auto iter = PII.begin(); iter != PII.end() ; ++iter)
	{
		DEBUG( *iter);
	}
	DEBUG("FILTERED");
	int count = 0;
//	map<int, int> real_id;
//	map<VertexId, int> old_vertices;
	for(auto iter = PIIFilter.begin(); iter != PIIFilter.end() ; ++iter)
	{
//		if ((old_vertices.find(g.EdgeEnd(iter->start)) != old_vertices.end()) && g.length(iter->start) > cfg::get().ds.IS) {
//			real_id.insert(make_pair(iter->id, iter->id));
//		} else
		{
			VertexId v = new_graph.AddVertex();
			newIds.AddVertexIntId(v, iter->id);
//			real_id.insert(make_pair(iter->id, iter->id));
//			old_vertices.insert(make_pair(g.EdgeEnd(iter->start), iter->id));
		}
	}
	DEBUG("PahtSetNumber is "<< PIIFilter.size());
	for(auto iter = PIIFilter.begin(); iter != PIIFilter.end() ; ++iter)
	{
		DEBUG( *iter);
	}

	map<EdgeId, double> weight_sums;
	for(auto iter = PIIFilter.begin(); iter != PIIFilter.end() ; ++iter)
	{
		PathSet<EdgeId> first = *iter;
		EdgeId old_first_edge = iter->start;
		if (weight_sums.find(old_first_edge) == weight_sums.end())
			weight_sums.insert (make_pair(old_first_edge, 0));
		else {
			INFO(weight_sums[old_first_edge] <<" " <<first.weight <<" "<< g.length(old_first_edge));
		}
		vector<PathSet<EdgeId>> extends;
		PI.FindExtension(PIIFilter,first, extends);
		if (extends.size() == 0){
			auto current_path = first.paths.begin();
			for(auto path_iter = current_path->begin(); path_iter != current_path->end(); ++path_iter) {
				old_first_edge = * path_iter;
				if (weight_sums.find(old_first_edge) == weight_sums.end())
					weight_sums.insert (make_pair(old_first_edge, 0));
				else {
					INFO(weight_sums[old_first_edge] <<" " <<first.weight)
				}
				weight_sums[old_first_edge] += first.weight;
			}
			old_first_edge = first.end;
			if (weight_sums.find(old_first_edge) == weight_sums.end())
				weight_sums.insert (make_pair(old_first_edge, 0));
			weight_sums[old_first_edge] += first.weight;
		} else {
			vector<PathSet<EdgeId>> extends;
			PI.FindExtension(PIIFilter,first, extends);
			for(size_t i = 0; i < extends.size(); i++)
			weight_sums[old_first_edge] += extends[i].weight;
		}
	}

	for(auto iter = PIIFilter.begin(); iter != PIIFilter.end() ; ++iter)
	{
		PathSet<EdgeId> first = *iter;
		vector<PathSet<EdgeId>> extends;
		PI.FindExtension(PIIFilter,first, extends);
		VertexId new_start = newIds.ReturnVertexId(iter->id);
		EdgeId old_first_edge = iter->start;
		if (first.weight / weight_sums[old_first_edge] < 0.9)
				DEBUG("low covered");
		DEBUG ("path-set numero " << first.id<< " has "<< extends.size()<<"extensions: ");
		for(size_t i = 0; i < extends.size(); i++) {
			DEBUG("to pathset "<< extends[i].id << " weight "<< extends[i].weight);
			VertexId new_end = newIds.ReturnVertexId(extends[i].id);
			DEBUG("adding edge from" << newIds.ReturnIntId(new_start) << " to " << newIds.ReturnIntId(new_end) << " of length " << g.length(old_first_edge) <<" and coverage "<< g.coverage(old_first_edge) << " * " << extends[i].weight / weight_sums[old_first_edge]);
			EdgeId eid = new_graph.AddEdge(new_start, new_end, g.EdgeNucls(old_first_edge));
			WrappedSetCoverage(new_graph, eid, (int) (g.coverage(old_first_edge) * g.length(old_first_edge) * extends[i].weight / weight_sums[old_first_edge]));
			DEBUG("count was "<< count);
//		    omnigraph::WriteSimple(new_graph, tot_labeler_after, cfg::get().output_dir  + ToString(count)+".dot", "no_repeat_graph");
			count ++ ;
		}
		if (extends.size() == 0){
			VertexId new_end = new_graph.AddVertex();
			DEBUG("adding edge from" << newIds.ReturnIntId(new_start) << " to " << newIds.ReturnIntId(new_end) << " of length " << g.length(old_first_edge) << " and coverage"<< g.coverage(old_first_edge) /*<< " * "  << first.weight / weight_sums[old_first_edge]*/);

			old_first_edge = first.start;
			EdgeId eid = new_graph.AddEdge(new_start, new_end, g.EdgeNucls(old_first_edge));
			WrappedSetCoverage(new_graph, eid, (int) (g.coverage(old_first_edge) * g.length(old_first_edge) /** first.weight / weight_sums[old_first_edge]*/));
			new_start = new_end;
			count++;

			if (first.paths.size() != 1 ) {
				DEBUG ("can not select between different tails:(");
			} else {
				DEBUG("singleton dead end!");
				auto current_path = first.paths.begin();
				DEBUG(current_path->size());
			//	DEBUG(first);


				for(auto path_iter = current_path->begin(); path_iter != current_path->end(); ++path_iter) {

					VertexId new_end = new_graph.AddVertex();
					old_first_edge = * path_iter;

					DEBUG("adding edge from" << newIds.ReturnIntId(new_start) << " to " << newIds.ReturnIntId(new_end) << " of length " << g.length(old_first_edge) << " and coverage"<< g.coverage(old_first_edge) /*<< " * "  << first.weight / weight_sums[old_first_edge]*/);

					EdgeId eid = new_graph.AddEdge(new_start, new_end, g.EdgeNucls(old_first_edge));
					WrappedSetCoverage(new_graph, eid, (int) (g.coverage(old_first_edge) * g.length(old_first_edge)));
					new_start = new_end;
					count++;

				}
				VertexId new_end = new_graph.AddVertex();
				old_first_edge = first.end;
				DEBUG("adding edge from" << newIds.ReturnIntId(new_start) << " to " << newIds.ReturnIntId(new_end) << " of length " << g.length(old_first_edge));


				EdgeId eid = new_graph.AddEdge(new_start, new_end, g.EdgeNucls(old_first_edge));
				WrappedSetCoverage(new_graph, eid, (int) (g.coverage(old_first_edge)  * g.length(old_first_edge) /** first.weight / weight_sums[old_first_edge] */));
				new_start = new_end;
				DEBUG("and tail of length "<< iter->length);
//				omnigraph::WriteSimple(new_graph, tot_labeler_after, cfg::get().output_dir  + ToString(count)+".dot", "no_repeat_graph");
				count++;
			}
		}
	}
	DEBUG(count);

}


};
}
