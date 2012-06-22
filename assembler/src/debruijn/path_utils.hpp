/*
 * path_utils.hpp
 *
 */

#pragma once

namespace debruijn_graph {

template<class Graph>
vector<typename Graph::EdgeId> GetCommonPathsEnd(
		const Graph& g,
		typename Graph::EdgeId& first_edge,
		typename Graph::EdgeId& second_edge, size_t min_dist,
		size_t max_dist) {
	PathStorageCallback<Graph> callback(g);
	PathProcessor<Graph> path_processor(g,
			min_dist - g.length(first_edge),
			max_dist - g.length(first_edge),
			g.EdgeEnd(first_edge), g.EdgeStart(second_edge),
			callback);
	path_processor.Process();
	auto paths = callback.paths();
	vector<typename Graph::EdgeId> result;
	if (paths.size() == 0) return result;
	if (paths.size() == 1) return paths[0];
	size_t j=0;
	while (j<paths[0].size()) {
		for (size_t i = 1;  i < paths.size(); ++i){
			if (j == paths[i].size()) {
				vector<typename Graph::EdgeId> result(paths[0].begin()+(paths[0].size() - j), paths[0].end());
				return result;
			} else {
				if (paths[0][paths[0].size()-1-j] != paths[i][paths[i].size()-1-j]) {
					vector<typename Graph::EdgeId> result(paths[0].begin()+(paths[0].size() - j), paths[0].end());
					return result;
				}
			}
		}
		j++;
	}
	return paths[0];

}



template<class Graph>
vector<vector<typename Graph::EdgeId> > GetAllPathsBetweenEdges(
		const Graph& g,
		typename Graph::EdgeId& first_edge,
		typename Graph::EdgeId& second_edge, size_t min_dist,
		size_t max_dist) {
	PathStorageCallback<Graph> callback(g);
	PathProcessor<Graph> path_processor(g,
			min_dist,
			max_dist, //0, *cfg::get().ds.IS - K + size_t(*cfg::get().ds.is_var),
			g.EdgeEnd(first_edge), g.EdgeStart(second_edge),
			callback);
	path_processor.Process();
	auto paths = callback.paths();
	return paths;
}

template<class graph_pack>
size_t GetAllPathsQuantity(const graph_pack& origin_gp,
		typename graph_pack::graph_t::EdgeId& first_edge,
		typename graph_pack::graph_t::EdgeId& second_edge, double dist) {
	PathStorageCallback<typename graph_pack::graph_t> callback(origin_gp.g);
	PathProcessor<typename graph_pack::graph_t> path_processor(
			origin_gp.g,
			dist - origin_gp.g.length(first_edge)
					- size_t(*cfg::get().ds.is_var),
			dist - origin_gp.g.length(first_edge)
					+ size_t(*cfg::get().ds.is_var),
			origin_gp.g.EdgeEnd(first_edge), origin_gp.g.EdgeStart(second_edge),
			callback);
	path_processor.Process();
	auto paths = callback.paths();
	TRACE(
			origin_gp.int_ids.ReturnIntId(first_edge) << " "
					<< origin_gp.int_ids.ReturnIntId(second_edge) << " "
					<< paths.size());
	return paths.size();
}


template<class graph_pack>
void GenerateMatePairStats(const graph_pack& origin_gp,
		PairedInfoIndex<typename graph_pack::graph_t>& clustered_index) {
	map < size_t, size_t > sizes;
	for (auto e_iter = origin_gp.g.SmartEdgeBegin(); !e_iter.IsEnd();
			++e_iter) {
		auto pi = clustered_index.GetEdgeInfo(*e_iter);
		for (auto i_iter = pi.begin(); i_iter != pi.end(); ++i_iter) {
			if (i_iter->d >= origin_gp.g.length(i_iter->first)) {
				size_t tmp = GetAllPathsQuantity(origin_gp, i_iter->first,
						i_iter->second, i_iter->d);
				if (sizes.find(tmp) == sizes.end())
					sizes.insert(make_pair(tmp, 0));
				sizes[tmp]++;
			}
		}
	}
	INFO("Pathset mate pair statistics:");
	for (auto s_iter = sizes.begin(); s_iter != sizes.end(); s_iter++) {
		INFO("- size: " << s_iter->first << "; pathsets: " << s_iter->second);
	}
}

}
