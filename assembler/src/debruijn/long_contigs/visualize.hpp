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
class PathsGraphLabeler : public GraphLabeler<Graph> {

protected:
	typedef GraphLabeler<Graph> super;
	typedef typename super::EdgeId EdgeId;
	typedef typename super::VertexId VertexId;
	Graph& g_;
	std::vector<BidirectionalPath>& paths_;
	std::map<EdgeId, std::string> labels_;

public:
	PathsGraphLabeler(Graph& g, std::vector<BidirectionalPath>& paths) : g_(g), paths_(paths){
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

void WriteGraphWithPathsSimple(const string& file_name, const string& graph_name, Graph& g, std::vector<BidirectionalPath>& paths,
		Path<Graph::EdgeId> path1, Path<Graph::EdgeId> path2) {

	std::fstream filestr;
	filestr.open(file_name.c_str(), std::fstream::out);

	gvis::DotGraphPrinter<Graph::VertexId> gp(graph_name, filestr);
	PathColorer<Graph> path_colorer(g, path1, path2);
	map<Graph::EdgeId, string> coloring = path_colorer.ColorPath();

	PathsGraphLabeler<Graph> gl(g, paths);
	ColoredGraphVisualizer<Graph> gv(g, gp, gl, coloring);
	AdapterGraphVisualizer<Graph> result_vis(g, gv);
	result_vis.Visualize();
	filestr.close();
}

} // namespace long_contigs

#endif /* VISUALIZE_HPP_ */
