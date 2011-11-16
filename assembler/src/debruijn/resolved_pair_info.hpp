#pragma once

#include "omni/edge_labels_handler.hpp"
#include "omni/paired_info.hpp"

namespace debruijn_graph {

template<class Graph>
class ResolvedGraphPairInfoCounter {
	typedef typename Graph::EdgeId EdgeId;

	const Graph& old_graph_;
	const PairedInfoIndex<Graph>& old_pair_info_;
	const Graph& new_graph_;
	const EdgeLabelHandler<Graph>& labels_;
public:
	ResolvedGraphPairInfoCounter(Graph& old_graph,
			PairedInfoIndex<Graph>& old_pair_info, Graph& new_graph,
			EdgeLabelHandler<Graph>& labels) :
			old_graph_(old_graph), old_pair_info_(old_pair_info), new_graph_(
					new_graph), labels_(labels) {
	}

	bool MapsUniquely(EdgeId old_edge) {
		VERIFY(
				labels_.edge_inclusions.find(old_edge)
						!= labels_.edge_inclusions.end());
//		VERIFY(labels_.edge_inclusions.find(old_edge)->second.size() > 0);
		if (!(labels_.edge_inclusions.find(old_edge)->second.size() > 0)){
			WARN("There are not exist copy for old graph edge "<< old_edge<<" "<<old_graph_.str(old_edge));
		}

		return labels_.edge_inclusions.find(old_edge)->second.size() == 1;
	}

	EdgeId UniqueMapping(EdgeId old_edge) {
		VERIFY(MapsUniquely(old_edge));
		return *(labels_.edge_inclusions.find(old_edge)->second.begin());
	}

	pair<EdgeId, size_t> OldEdgePositionInNewGraph(EdgeId old_edge) {
		VERIFY(MapsUniquely(old_edge));
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

	void ProcessInfos(PairedInfoIndex<Graph>& new_pair_info, const vector<PairInfo<EdgeId>>& infos) {
		EdgeId first_old_edge = infos.begin()->first;
		EdgeId second_old_edge = infos.begin()->second;
		if (MapsUniquely(first_old_edge) && MapsUniquely(second_old_edge)) {
			pair<EdgeId, size_t> new_pos_of_first = OldEdgePositionInNewGraph(first_old_edge);
			pair<EdgeId, size_t> new_pos_of_second = OldEdgePositionInNewGraph(second_old_edge);
			for (auto it = infos.begin(); it != infos.end(); ++it) {
				new_pair_info.AddPairInfo(
						PairInfo<EdgeId>(new_pos_of_first.first, new_pos_of_second.first,
								0. + it->d + new_pos_of_first.second - new_pos_of_second.second,
								it->weight, it->variance), false);
			}
		}
	}

	void FillResolvedGraphPairedInfo(PairedInfoIndex<Graph>& new_pair_info) {
		for (auto it = old_pair_info_.begin(); it != old_pair_info_.end(); ++it) {
			ProcessInfos(new_pair_info, *it);
		}
	}
};

}
