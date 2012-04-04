//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * visualize.hpp
 *
 *  Created on: Aug 4, 2011
 *      Author: andrey
 */

#ifndef VISUALIZE_HPP_
#define VISUALIZE_HPP_

#include "lc_common.hpp"

namespace long_contigs {

using namespace debruijn_graph;

template<class Graph>
class PathsGraphLabeler : public AbstractGraphLabeler<Graph> {
	typedef AbstractGraphLabeler<Graph> base;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	const std::vector<BidirectionalPath>& paths_;
	map<EdgeId, std::string> labels_;

public:
	PathsGraphLabeler(const Graph& g, const std::vector<BidirectionalPath>& paths) : base(g), paths_(paths) {
//		for (auto iter = g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
//			labels_[*iter] = "";
//		}

		for(size_t i = 0; i < paths.size(); ++i) {
			const BidirectionalPath& path = paths[i];

			for (size_t idx = 0; idx < path.size(); ++idx) {
				if (labels_[path[idx]].size() > 0) {
					labels_[path[idx]] += ", ";
				}
				labels_[path[idx]] += "(" + ToString(paths[i].uid) + " : " + ToString(idx) + ")";
			}
		}
	}

	virtual std::string label(VertexId vertexId) const {
		return "";
	}

	virtual std::string label(EdgeId edgeId) const {
		if(labels_.count(edgeId) == 0)
			return "";
		else
			return labels_.find(edgeId)->second;
	}
};

void WritePathLocality(const conj_graph_pack& gp, const GraphLabeler<Graph>& labeler,
		const string& folder, const BidirectionalPath& path, size_t edge_split_length
		, const Path<EdgeId>& color1, const Path<EdgeId>& color2) {
	WriteComponentsAlongPath(gp.g, labeler, folder + ToString(path.uid) + ".dot"
			, edge_split_length
			, as_trivial_mapping_path(gp.g, as_simple_path(gp.g, path))
			, color1, color2, true);
}

void WritePathLocalities(const conj_graph_pack& gp, const string& folder, const std::vector<BidirectionalPath>& paths) {
	auto path1 = FindGenomeMappingPath<K> (gp.genome, gp.g, gp.index, gp.kmer_mapper);
	auto path2 = FindGenomeMappingPath<K> (!gp.genome, gp.g, gp.index, gp.kmer_mapper);

	for (auto it = paths.begin(); it != paths.end(); ++it) {
		if (it->size() > 1) {
			/////todo code duplication
			StrGraphLabeler<Graph> str_labeler(gp.g);
			PathsGraphLabeler<Graph> path_labeler(gp.g, vector<BidirectionalPath>{*it});
			EdgePosGraphLabeler<Graph> pos_labeler(gp.g, gp.edge_pos);

			CompositeLabeler<Graph> composite_labeler(str_labeler, path_labeler, pos_labeler);
			/////todo code duplication

			string path_folder = folder + ToString(it->uid) + "/";
			make_dir(path_folder);
			WritePathLocality(gp, composite_labeler, path_folder, *it, 1000, path1.simple_path(), path2.simple_path());
		}
	}
}

void WriteGraphWithPathsSimple(const conj_graph_pack& gp, const string& file_name, const string& graph_name, const std::vector<BidirectionalPath>& paths) {
	std::fstream filestr;
	filestr.open(file_name.c_str(), std::fstream::out);


	StrGraphLabeler<Graph> str_labeler(gp.g);
	PathsGraphLabeler<Graph> path_labeler(gp.g, paths);
	EdgePosGraphLabeler<Graph> pos_labeler(gp.g, gp.edge_pos);

	CompositeLabeler<Graph> composite_labeler(str_labeler, path_labeler, pos_labeler);

	auto_ptr<GraphColorer<Graph>> colorer(DefaultColorer(gp.g, FindGenomePath<K> (gp.genome, gp.g, gp.index)
			, FindGenomePath<K> (!gp.genome, gp.g, gp.index)));

	DotGraphPrinter<Graph> printer(gp.g, composite_labeler, *colorer, graph_name, filestr);

	ColoredGraphVisualizer<Graph> gv(gp.g, printer);
	AdapterGraphVisualizer<Graph> result_vis(gp.g, gv);
	result_vis.Visualize();
	filestr.close();
}

} // namespace long_contigs

#endif /* VISUALIZE_HPP_ */
