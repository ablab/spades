//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "omni/edge_labels_handler.hpp"
#include "de/paired_info.hpp"

namespace debruijn_graph {

template<class Graph>
class ResolvedGraphPairInfoCounter {
	typedef typename Graph::EdgeId EdgeId;
  typedef set<Point> Histogram;

	const Graph& old_graph_;
	const PairedInfoIndexT<Graph>& old_pair_info_;
	const Graph& new_graph_;
	const EdgeLabelHandler<Graph>& labels_;

	size_t no_current_copies;
	set<int> no_current_copies_ids;

public:
	ResolvedGraphPairInfoCounter(Graph& old_graph,
			const PairedInfoIndexT<Graph>& old_pair_info, Graph& new_graph,
			EdgeLabelHandler<Graph>& labels) :
			old_graph_(old_graph), old_pair_info_(old_pair_info), new_graph_(
					new_graph), labels_(labels) {
	}

	bool MapsUniquely(EdgeId old_edge) {
		VERIFY_MSG(labels_.edge_inclusions.find(old_edge) != labels_.edge_inclusions.end(), 
               "old_edge " << old_graph_.int_id(old_edge));
//		VERIFY(labels_.edge_inclusions.find(old_edge)->second.size() > 0);
		if (!(labels_.edge_inclusions.find(old_edge)->second.size() > 0)){
			DEBUG("There are no current copy for old graph edge " << old_edge << " " << old_graph_.str(old_edge));
			no_current_copies_ids.insert((int) old_graph_.int_id(old_edge));
			no_current_copies++;
		}
		return labels_.edge_inclusions.find(old_edge)->second.size() == 1;
	}

	EdgeId UniqueMapping(EdgeId old_edge) {
		bool flag = MapsUniquely(old_edge);
		VERIFY(flag);
		return *(labels_.edge_inclusions.find(old_edge)->second.begin());
	}

	pair<EdgeId, size_t> OldEdgePositionInNewGraph(EdgeId old_edge) {
        bool flag = MapsUniquely(old_edge);
        VERIFY(flag);
		EdgeId new_edge = UniqueMapping(old_edge);
		VERIFY(labels_.edge_labels.find(new_edge) != labels_.edge_labels.end());
		const vector<EdgeId>& old_edges = labels_.edge_labels.find(new_edge)->second;
		VERIFY(std::find(old_edges.begin(), old_edges.end(), old_edge) != old_edges.end());
		size_t offset = 0;
		for (auto it = old_edges.begin(); it != old_edges.end(); ++it) {
			if (*it == old_edge)
			break;
			offset += old_graph_.length(*it);
		}
		return make_pair(new_edge, offset);
	}

	void ProcessInfos(PairedInfoIndexT<Graph>& new_pair_info, 
                    EdgeId e1, 
                    EdgeId e2, 
                    const Histogram& infos) 
  {
		if (MapsUniquely(e1) && MapsUniquely(e2)) {
			pair<EdgeId, size_t> new_pos_of_first = OldEdgePositionInNewGraph(e1);
			pair<EdgeId, size_t> new_pos_of_second = OldEdgePositionInNewGraph(e2);
			for (auto it = infos.begin(); it != infos.end(); ++it) {
				new_pair_info.AddPairInfo(new_pos_of_first.first, new_pos_of_second.first,
								0. + it->d + (double) new_pos_of_first.second - (double) new_pos_of_second.second,
								it->weight, it->var, false);
			}
		}
	}

	void FillResolvedGraphPairedInfo(PairedInfoIndexT<Graph>& new_pair_info) {
		no_current_copies = 0;
		no_current_copies_ids.clear();
		for (auto it = old_pair_info_.begin(); it != old_pair_info_.end(); ++it) {
			ProcessInfos(new_pair_info, it.first(), it.second(), *it);
		}
		if (no_current_copies) {
			DEBUG("There were no current copies for " << no_current_copies_ids.size() << " old graph edges: "<<ToString(no_current_copies_ids));
		}
	}

 private:
  DECL_LOGGER("ResolvedGraphPairInfoCounter")

};



}
