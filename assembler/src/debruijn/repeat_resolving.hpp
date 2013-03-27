//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "repeat_resolver.hpp"
#include "graphio.hpp"
#include "omni/loop_resolver.hpp"
#include <sys/types.h>
#include <sys/stat.h>

namespace debruijn_graph {

template<class Graph>
void ResolveRepeats(Graph &g, IdTrackHandler<Graph> &old_ids,
		const PairedInfoIndexT<Graph> &info,
		EdgesPositionHandler<Graph> &edges_pos, Graph &new_graph,
		IdTrackHandler<Graph> &new_ids,
		EdgesPositionHandler<Graph> &edges_pos_new, const string& output_folder,
		EdgeLabelHandler<Graph> &LabelsAfter, bool developer_mode) {
	INFO("SUBSTAGE == Resolving primitive repeats");
	for (auto e_iter = g.SmartEdgeBegin(); !e_iter.IsEnd(); ++e_iter) {
		{
			if ((g.EdgeStart(*e_iter) == g.EdgeEnd(*e_iter))) {
				if ((g.length(*e_iter) > 1)) {
					size_t half_len = g.length(*e_iter) / 2;
					g.SplitEdge(*e_iter, half_len);
				} else {
					WARN("Loop of length 1 detected");
				}
			}
		}
	}
	DeletedVertexHandler<Graph> tmp_deleted_handler(new_graph);
	TRACE("deleted handler created");
	RepeatResolver<Graph> repeat_resolver(g, old_ids, info, edges_pos,
			new_graph, new_ids, edges_pos_new, tmp_deleted_handler, LabelsAfter,
			developer_mode);
	make_dir(output_folder);
	auto edge_labels = repeat_resolver.GetEdgeLabels();
	LabelsAfter.FillLabels(edge_labels);

	repeat_resolver.ResolveRepeats(output_folder);
	DEBUG("Primitive repeats resolved");
}

}
