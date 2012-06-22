/*
 * split_path_constructor.hpp
 *
 *  Created on: Jun 14, 2012
 *      Author: avsirotkin
 */

#pragma once

#include "standard.hpp"
#include "logger/logger.hpp"
#include "path_utils.hpp"


namespace debruijn_graph {


template<class Graph>
class PathInfoClass {
public:
	typedef typename Graph::EdgeId EdgeId;
	typedef omnigraph::PairInfo<EdgeId> PairInfo;

	EdgeId base_edge;
	vector<PairInfo> path;
	PathInfoClass(): base_edge(NULL) {};
	PathInfoClass(const EdgeId Edge): base_edge(Edge) {};
	pair<EdgeId, double> operator[](const size_t i) const {
		if (i == 0) {
			return(make_pair(base_edge, 0.0));
		}
		VERIFY(i < path.size() + 1);
		return (make_pair(path[i-1].second, path[i-1].d));
	}
	size_t size() {
		return path.size() + 1;
	}
	void push_back(const PairInfo& pi) {
		path.push_back(pi);
	}
	typename vector<PairInfo>::iterator begin() {
		return path.begin();
	}
	typename vector<PairInfo>::iterator end() {
		return path.end();
	}
	std::string PrintPath(const Graph& graph) {
		std::ostringstream ss;
		ss<<" "<<graph.int_id(base_edge)<<": ";
		for (size_t j=0; j < path.size(); j++){
			ss<<"("<<graph.int_id(path[j].second)<<", "<<path[j].d<<"), ";
		}
		return ss.str();
	}

};





template<class Graph>
class SplitPathConstructor {
	const Graph &graph_;
	map<EdgeId, vector<PathInfoClass<Graph> > > split_pathes;
public:
	SplitPathConstructor(const Graph &graph): graph_(graph) {};
	vector<PathInfoClass<Graph>> ConvertEdgePairInfoToSplitPathes(vector<PairInfo<typename Graph::EdgeId>> &p_infos){
		EdgeId cur_edge;
		vector<PathInfoClass<Graph>> result;
		if (p_infos.size() == 0) return result;
		else {
			cur_edge = p_infos[0].first;
			for (size_t i = 1; i < p_infos.size(); ++i)
				if (cur_edge != p_infos[i].first){
					WARN("Can not convert pair infos from different edges into split paths");
					VERIFY(false);
				}
			set<size_t> pair_info_used;
			for (size_t i = p_infos.size(); i > 0; --i){
				if (p_infos[i-1].d < 0.01) continue;
				if (pair_info_used.find(i-1) != pair_info_used.end()) continue;
				DEBUG("SPC: pi "<<graph_.int_id(cur_edge)<<" "<<graph_.int_id(p_infos[i-1].second)<<" "<<p_infos[i-1].d<<" "<<p_infos[i-1].variance);
				auto common_part = GetCommonPathsEnd(graph_, cur_edge, p_infos[i-1].second, p_infos[i-1].d - p_infos[i-1].variance, p_infos[i-1].d + p_infos[i-1].variance);
				DEBUG("Found common part of size "<<common_part.size());
				PathInfoClass<Graph> sub_res(cur_edge);
				if (common_part.size() > 0) {
					size_t total_length = 0;
					for (size_t j = 0; j < common_part.size(); j++){
						total_length += graph_.length(common_part[j]);
					}
					DEBUG(ToString(common_part));
					for (size_t j = 0; j < common_part.size(); j++){
						PairInfo<typename Graph::EdgeId> cur_pi(cur_edge, common_part[j], p_infos[i-1].d - total_length, p_infos[i-1].weight, p_infos[i-1].variance);

						sub_res.push_back(cur_pi);
						total_length -= graph_.length(common_part[j]);
						for (size_t ind = 0; ind < i-1; ind++){
							if (cur_pi.first == p_infos[ind].first
									&&  cur_pi.second == p_infos[ind].second
									&&  cur_pi.d == p_infos[ind].d) {
								pair_info_used.insert(ind);
							}

						}

					}
				}
				sub_res.push_back(p_infos[i-1]);
				result.push_back(sub_res);
				DEBUG(sub_res.PrintPath(graph_));
			}
		}
		return result;
	}
};


}

