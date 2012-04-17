//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once
#include "matepair_transformer.hpp"
#include "path_set.hpp"
#include "new_debruijn.hpp"
#include "graph_pack.hpp"
#include "utils.hpp"
#include "omni/omni_utils.hpp"
#include "verify.hpp"
#include "omni/omni_tools.hpp"
#include "omni/omnigraph.hpp"
#include <unordered_map>
#include "omni/edges_position_handler.hpp"
#include "omni/total_labeler.hpp"
#include "path_set_stats.hpp"
#include "path_set_tools.hpp"
namespace debruijn_graph{
template <class graph_pack>
class PathSetGraphConstructor {

typedef typename graph_pack::graph_t::EdgeId EdgeId;
typedef typename graph_pack::graph_t::VertexId VertexId;
typedef typename graph_pack::graph_t Graph;
typedef vector<EdgeId > Path;
graph_pack& gp;
const PairedInfoIndex<Graph>& pair_info_;
graph_pack& new_gp;
RestrictedMap<EdgeId, EdgeId> new_to_old;
RestrictedMap<EdgeId, EdgeId> old_to_new;

public:
	PathSetGraphConstructor(graph_pack& gp, PairedInfoIndex<Graph>& clustered_index, graph_pack& new_gp): gp(gp), pair_info_(clustered_index), new_gp(new_gp)
	{
		new_to_old.clear();
		old_to_new.clear();

	}

	void ConstructOverlap() {
		PathSetIndexData<EdgeId> PII ;
		PathSetIndexData<EdgeId> PIIFilter ;

		MatePairTransformer<graph_pack> transformer(gp, pair_info_);
		transformer.Transform(PII);
		PathSetIndex<graph_pack> PI(PII, gp);
		//	PI.RemovePrefixes(PIIFilter_tmp);
		PI.Process(PIIFilter);


		for(auto iter = PIIFilter.begin(); iter != PIIFilter.end() ; ++iter)
		{
			INFO(str(*iter, gp));
			//		DEBUG(tst());
		}
		DEBUG("FILTERED");
//		int count = 0;
		map<PathSet<EdgeId>, vector<PathSet<EdgeId> > > extentionMap, backwardMap;

		for(auto iter = PIIFilter.begin(); iter != PIIFilter.end() ; ++iter) {
			vector<PathSet<EdgeId> > extends;
			PathSet<EdgeId> first = *iter;

			PI.FindExtension(PIIFilter, first, extends);
			extentionMap[*iter] = extends;
			for(auto inside_iter = extends.begin(); inside_iter != extends.end(); ++inside_iter) {
				if (backwardMap.find(*inside_iter) == backwardMap.end()) {
					vector<PathSet<EdgeId> > aux;
					backwardMap[*inside_iter] = aux;
				}
				backwardMap[*inside_iter].push_back(*iter);
			}

			//		DEBUG(tst());
		}
		for(auto iter = PIIFilter.begin(); iter != PIIFilter.end() ; ++iter) {
			DEBUG ("id: " << iter->id<< " fwd " <<extentionMap[*iter].size());
			for (auto tmp_iter = extentionMap[*iter].begin(); tmp_iter != extentionMap[*iter].end(); tmp_iter++) {
				DEBUG(tmp_iter->id);
			}
			DEBUG( " bwd "<< backwardMap[*iter].size());
			for (auto tmp_iter = backwardMap[*iter].begin(); tmp_iter != backwardMap[*iter].end(); tmp_iter++) {
				DEBUG(tmp_iter->id);
			}
		}
		for(auto iter = PIIFilter.begin(); iter != PIIFilter.end() ; ++iter) {
			bool strange_start = ((extentionMap[*iter].size() == 1) && (backwardMap[*iter].size() == 1) && (extentionMap[backwardMap[*iter][0]].size() > 1));
			if (((extentionMap[*iter].size() <= 1) && (backwardMap[*iter].size() != 1)) || strange_start) {
				if(strange_start) {
					DEBUG("strange start from pathset id " << iter-> id <<" edge "<< gp.int_ids.ReturnIntId(iter->start));
				}

				VertexId v = new_gp.g.AddVertex();
				new_gp.int_ids.AddVertexIntId(v, iter->id);
				DEBUG("working on path starting from " << iter->id);
				auto tmp_iter = *iter;
				while (extentionMap[tmp_iter].size() == 1) {
					VertexId end = new_gp.g.AddVertex();
					new_gp.int_ids.AddVertexIntId(end, extentionMap[tmp_iter][0].id);
					AddEdgeWithAllHandlers(v, end , tmp_iter.start);
					v = end;
					tmp_iter = extentionMap[tmp_iter][0];
					if (backwardMap[tmp_iter].size() != 1)
						break;
				}
				if (extentionMap[tmp_iter].size() == 0) {
					if (tmp_iter.paths.size() != 1) {
						WARN("Non unique tail, removed");
					} else {
						DEBUG("tail adding..");
						auto current_path = tmp_iter.paths.begin();
						DEBUG("..of length: "<<current_path->size());
						VertexId end = new_gp.g.AddVertex();
						new_gp.int_ids.AddVertexIntId(end, - new_gp.int_ids.ReturnIntId(end));
						EdgeId old_edge = tmp_iter.start;
						AddEdgeWithAllHandlers(v, end , old_edge);
						v = end;
						for(auto path_iter = current_path->begin(); path_iter != current_path->end(); ++path_iter) {
							end = new_gp.g.AddVertex();
							new_gp.int_ids.AddVertexIntId(end, - new_gp.int_ids.ReturnIntId(end));
							old_edge = *path_iter;
							AddEdgeWithAllHandlers(v, end , old_edge);
							v = end;
						}
						end = new_gp.g.AddVertex();
						new_gp.int_ids.AddVertexIntId(end, - new_gp.int_ids.ReturnIntId(end));
						old_edge = tmp_iter.end;
						AddEdgeWithAllHandlers(v, end , old_edge);
						v = end;

						DEBUG("tail added");
					}
				} else {
					DEBUG("came into many-one vertex");
				}
			}
		}


		INFO("Adding many-outgoing edges");
		for(auto iter = PIIFilter.begin(); iter != PIIFilter.end() ; ++iter)
		{
			//if edge was added as a tail because of gap in coverage
			if ((extentionMap[*iter].size() > 1)) {
				INFO ("id: " << iter->id<< " fwd " <<extentionMap[*iter].size() << " bwd "<< backwardMap[*iter].size() << "  extentions: ");
				if (old_to_new.find(iter->start) != old_to_new.end()) {
					INFO("Edge" << gp.int_ids.ReturnIntId(iter->start) <<" already added!(as a tail?)");
				} else {
					for (auto tmp_iter = extentionMap[*iter].begin(); tmp_iter < extentionMap[*iter].end(); ++tmp_iter) {
						INFO(tmp_iter->id);
					}
					VertexId v = new_gp.g.AddVertex();
					new_gp.int_ids.AddVertexIntId(v, iter->id);
					VertexId end = new_gp.g.AddVertex();
					new_gp.int_ids.AddVertexIntId(end, - 10000 -new_gp.int_ids.ReturnIntId(end));
					AddEdgeWithAllHandlers(v, end , iter->start);
				}
			}
		}
		//	map<VertexId, int> long_start_verticesmaxId;


		INFO("adding isolated edges");

		for(auto iter = gp.g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
			VertexId start = gp.g.EdgeStart(*iter);
			VertexId end = gp.g.EdgeEnd(*iter);
			int flag = 0;
			TRACE (gp.g.CheckUniqueOutgoingEdge(start)<<" "<<  gp.g.IsDeadStart(start) <<" "<< gp.g.CheckUniqueIncomingEdge(end) <<" "<<gp.g.IsDeadEnd(end));
			if (old_to_new.find(*iter) == old_to_new.end())
				flag = 1;
			if ((gp.g.CheckUniqueOutgoingEdge(start) && gp.g.IsDeadStart(start) && gp.g.CheckUniqueIncomingEdge(end) && gp.g.IsDeadEnd(end)) || flag ) {
				if (flag) {
					WARN(" adding non isolated, non added edge: " << gp.int_ids.ReturnIntId(*iter));
				}
				VertexId new_start = new_gp.g.AddVertex();
				VertexId new_end = new_gp.g.AddVertex();
				new_gp.int_ids.AddVertexIntId(new_start, - 100000 -new_gp.int_ids.ReturnIntId(new_start));
				new_gp.int_ids.AddVertexIntId(new_end, - 100000 -new_gp.int_ids.ReturnIntId(new_end));

				EdgeId old_first_edge = *iter;
				AddEdgeWithAllHandlers(new_start, new_end, old_first_edge);
			}
		}

	}
	void Construct() {
		PathSetIndexData<EdgeId> PII ;
		PathSetIndexData<EdgeId> PIIFilter ;

		MatePairTransformer<graph_pack> transformer(gp, pair_info_);
		transformer.Transform(PII);
		PathSetIndex<graph_pack> PI(PII, gp);
		//	PI.RemovePrefixes(PIIFilter_tmp);
		PI.Process(PIIFilter);


		for(auto iter = PII.begin(); iter != PII.end() ; ++iter)
		{
			DEBUG(str(*iter, gp));
		}
		DEBUG("FILTERED");
		int count = 0;
		for(auto iter = PIIFilter.begin(); iter != PIIFilter.end() ; ++iter)
		{
			VertexId v = new_gp.g.AddVertex();
			new_gp.int_ids.AddVertexIntId(v, iter->id);
		}
		map<EdgeId, int> long_starts;
		map<int, int> real_ids;
		for(auto iter = PIIFilter.begin(); iter != PIIFilter.end() ; ++iter)
		{
			DEBUG("watching on id "<<  iter->id);
			if (long_starts.find(iter->start) == long_starts.end()) {
				real_ids.insert(make_pair(iter->id, iter->id));
				long_starts.insert(make_pair(iter->start, iter->id));
			} else {
				bool was_glued = false;
				DEBUG("interesting");
				for(auto iter2 = PIIFilter.begin(); iter2 != iter; ++iter2 ){
					if (iter2 == iter)
						break;
					if (NeedToGlue(*iter, *iter2)) {
						DEBUG("gluing");
						real_ids.insert(make_pair(iter->id, real_ids[iter2->id]));
						new_gp.g.DeleteVertex(new_gp.int_ids.ReturnVertexId(iter->id));

						was_glued = true;
						break;
					}
				}
				if(!was_glued) {
					real_ids.insert(make_pair(iter->id, iter->id));
				}
			}
		}

		DEBUG("PahtSetNumber is "<< PIIFilter.size());
		for(auto iter = PIIFilter.begin(); iter != PIIFilter.end() ; ++iter)
		{
			DEBUG(str(*iter, gp));
		}
		map<EdgeId, double> weight_sums;
		for(auto iter = PIIFilter.begin(); iter != PIIFilter.end() ; ++iter)
		{
			PathSet<EdgeId> first = *iter;
			EdgeId old_first_edge = iter->start;
			if (weight_sums.find(old_first_edge) == weight_sums.end())
				weight_sums.insert (make_pair(old_first_edge, 0));
			else {
				INFO(weight_sums[old_first_edge] <<" " <<first.weight <<" "<< gp.g.length(old_first_edge));
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
			VertexId new_start = new_gp.int_ids.ReturnVertexId(real_ids[iter->id]);
			EdgeId old_first_edge = iter->start;
			if (first.weight / weight_sums[old_first_edge] < 0.9)
				DEBUG("low covered");
			DEBUG ("path-set numero " << first.id<< " has "<< extends.size()<<"extensions: ");
			for(size_t i = 0; i < extends.size(); i++) {
				DEBUG("to pathset "<< extends[i].id << " and vertex "<<real_ids[extends[i].id] << " weight "<< extends[i].weight);
				VertexId new_end = new_gp.int_ids.ReturnVertexId(real_ids[extends[i].id]);
				VERIFY(new_end != VertexId(NULL));
				VERIFY(real_ids[extends[i].id] == new_gp.int_ids.ReturnIntId(new_end));
//				DEBUG("adding edge from" << new_gp.int_ids.ReturnIntId(new_start) << " to " << new_gp.int_ids.ReturnIntId(new_end) << " of length " << gp.g.length(old_first_edge) <<" and coverage "<< gp.g.coverage(old_first_edge) << " * " << extends[i].weight / weight_sums[old_first_edge]);
				bool flag = true;
				vector<EdgeId> out_e = new_gp.g.OutgoingEdges(new_start);
				for (auto eid = out_e.begin(); eid != out_e.end(); ++eid){
					if (new_to_old[*eid] == old_first_edge && new_gp.g.EdgeEnd(*eid)== new_end) {
						//TODO: what's about coverage?
						flag = false;
						break;
					}
				}
				if ( (!flag)/* && real_ids[first.id] != first.id*/) {
					DEBUG("ignoring clone to pathset " << extends[i].id << " and vertex " << real_ids[extends[i].id]);
				} else {
					AddEdgeWithAllHandlers(new_start, new_end, old_first_edge);
					count ++ ;
				}
			}
			if (extends.size() == 0){
				VertexId new_end = new_gp.g.AddVertex();
				DEBUG("adding edge from" << new_gp.int_ids.ReturnIntId(new_start) << " to " << new_gp.int_ids.ReturnIntId(new_end) << " of length " << gp.g.length(old_first_edge) << " and coverage"<< gp.g.coverage(old_first_edge) /*<< " * "  << first.weight / weight_sums[old_first_edge]*/);

				old_first_edge = first.start;
				AddEdgeWithAllHandlers(new_start, new_end, old_first_edge);

				new_start = new_end;
				count++;

				if (first.paths.size() != 1 ) {
					DEBUG ("can not select between different tails:(");
				} else {
					DEBUG("singleton dead end!");
					auto current_path = first.paths.begin();
					DEBUG(current_path->size());
					for(auto path_iter = current_path->begin(); path_iter != current_path->end(); ++path_iter) {
						VertexId new_end = new_gp.g.AddVertex();
						old_first_edge = * path_iter;
						AddEdgeWithAllHandlers(new_start, new_end, old_first_edge);
						count++;
					}
					VertexId new_end = new_gp.g.AddVertex();
					old_first_edge = first.end;
					AddEdgeWithAllHandlers(new_start, new_end, old_first_edge);

//					EdgeId eid = new_gp.g.AddEdge(new_start, new_end, gp.g.EdgeNucls(old_first_edge));
//					WrappedSetCoverage(new_gp.g, eid, (int) (gp.g.coverage(old_first_edge)  * gp.g.length(old_first_edge) /** first.weight / weight_sums[old_first_edge] */));
//					new_gp.edge_pos.AddEdgePosition(eid, gp.edge_pos.edges_positions().find(old_first_edge)->second);
//					new_to_old.insert(make_pair(eid, old_first_edge));
					new_start = new_end;
					DEBUG("and tail of length "<< iter->length);
					//				omnigraph::WriteSimple(new_gp.g, tot_labeler_after, cfg::get().output_dir  + ToString(count)+".dot", "no_repeat_graph");
					count++;
				}
			}
		}
		DEBUG(count);
		INFO("adding isolated edges");

		for(auto iter = gp.g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
			VertexId start = gp.g.EdgeStart(*iter);
			VertexId end = gp.g.EdgeEnd(*iter);
			TRACE (gp.g.CheckUniqueOutgoingEdge(start)<<" "<<  gp.g.IsDeadStart(start) <<" "<< gp.g.CheckUniqueIncomingEdge(end) <<" "<<gp.g.IsDeadEnd(end));
			if (gp.g.CheckUniqueOutgoingEdge(start) && gp.g.IsDeadStart(start) && gp.g.CheckUniqueIncomingEdge(end) && gp.g.IsDeadEnd(end) ) {
				VertexId new_start = new_gp.g.AddVertex();
				VertexId new_end = new_gp.g.AddVertex();
				EdgeId old_first_edge = *iter;
				AddEdgeWithAllHandlers(new_start, new_end, old_first_edge);
			}
		}
	}
	bool NeedToGlue(const PathSet<EdgeId>& e1,const PathSet<EdgeId>& e2){
		if (e1.start != e2.start) {
			return false;
		} else {
			bool glue = false;
			DEBUG("checking " << str(e1, gp) <<" and " <<str(e2, gp));
			for(auto iter1 = e1.paths.begin(); iter1 != e1.paths.end(); ++iter1) {
				for(auto iter2 = e2.paths.begin(); iter2 != e2.paths.end(); ++iter2) {
					size_t sum_len = gp.g.length(e1.start);
					size_t it = 0;
					if (iter2->size() != 0) {
						for(size_t iter = 0; iter!= iter1->size(); ++iter) {
							if ((*iter1)[iter] == (*iter2)[it]) {
								sum_len += gp.g.length((*iter1)[iter]);
								if (it != iter2->size()) {
									++it;
								} else {
									break;
								}
							} else {
								break;
							}
						}
					}
					DEBUG("sum len"<< sum_len << " " << cfg::get().ds.IS);
					if (sum_len > *cfg::get().ds.IS - cfg::get().rr.near_vertex)
						glue = true;
				}

			}
			DEBUG ("result is" << glue);
			return glue;
		}
	}


	RestrictedMap<EdgeId, EdgeId> GetEdgeLabels(){
		return new_to_old;
	}

	EdgeId AddEdgeWithAllHandlers(VertexId new_start, VertexId new_end, EdgeId old_first_edge){
		DEBUG("adding edge from" << new_gp.int_ids.ReturnIntId(new_start) << " to " << new_gp.int_ids.ReturnIntId(new_end) << " of length " << gp.g.length(old_first_edge));


		EdgeId eid = new_gp.g.AddEdge(new_start, new_end, gp.g.EdgeNucls(old_first_edge));
		WrappedSetCoverage(new_gp.g, eid, (int) (gp.g.coverage(old_first_edge)  * gp.g.length(old_first_edge) /** first.weight / weight_sums[old_first_edge] */));
		new_gp.edge_pos.AddEdgePosition(eid, gp.edge_pos.edges_positions().find(old_first_edge)->second);
		new_to_old.insert(make_pair(eid, old_first_edge));
		old_to_new.insert(make_pair( old_first_edge, eid));
		new_start = new_end;
		DEBUG("edge added");
		//DEBUG("and tail of length "<< iter->length);
		return eid;
	}


};
}
