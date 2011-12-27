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

	std::vector<BidirectionalPath>& paths_;
	std::map<EdgeId, std::string> labels_;

public:
	PathsGraphLabeler(const Graph& g, std::vector<BidirectionalPath>& paths) : base(g), paths_(paths){
		for(size_t i = 0; i < paths.size(); ++i) {
			BidirectionalPath& path = paths[i];

			for (auto iter = g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
				size_t edgeLen = g.length(*iter);
				std::string label = "{" + ToString(edgeLen) + "}";
				labels_.insert(std::make_pair(*iter, label));
			}

			for(size_t edge = 0; edge < path.size(); ++edge) {
				labels_[path[edge]] += "," + ToString(i) + "(" + ToString(edge) + ")";
			}
		}
	}

	virtual std::string label(VertexId vertexId) const {
		return "";
	}

	virtual std::string label(EdgeId edgeId) const {
		auto label = labels_.find(edgeId);
		return label == labels_.end() ? "" : label->second;
	}
};

void WriteGraphWithPathsSimple(const conj_graph_pack& gp, const string& file_name, const string& graph_name, std::vector<BidirectionalPath>& paths) {

	Path<Graph::EdgeId> path1 = FindGenomePath<K> (gp.genome, gp.g, gp.index);
	Path<Graph::EdgeId> path2 = FindGenomePath<K> (!gp.genome, gp.g, gp.index);

	std::fstream filestr;
	filestr.open(file_name.c_str(), std::fstream::out);

	gvis::DotGraphPrinter<Graph::VertexId> printer(graph_name, filestr);
	PathColorer<Graph> path_colorer(gp.g, path1, path2);

	map<Graph::EdgeId, string> coloring = path_colorer.ColorPath();

	StrGraphLabeler<Graph> str_labeler(gp.g);
	PathsGraphLabeler<Graph> path_labeler(gp.g, paths);
	EdgePosGraphLabeler<Graph> pos_labeler(gp.g, gp.edge_pos);

	CompositeLabeler<Graph> composite_labeler = {&str_labeler, &path_labeler, &pos_labeler};

	ColoredGraphVisualizer<Graph> gv(gp.g, printer, composite_labeler, coloring);
	AdapterGraphVisualizer<Graph> result_vis(gp.g, gv);
	result_vis.Visualize();
	filestr.close();
}

} // namespace long_contigs

#endif /* VISUALIZE_HPP_ */
