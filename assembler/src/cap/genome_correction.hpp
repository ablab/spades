#pragma once

#include "utils.hpp"
#include "coloring.hpp"
//#include <map>;

namespace cap {

template<class T>
class bag {
	/*std::*/
	map<T, size_t> data_;
public:
	typedef typename map<T, size_t>::const_iterator const_iterator;

	void put(const T& t, size_t mult) {
		VERIFY(mult > 0);
		data_[t] += mult;
	}

	void put(const T& t) {
		put(t, 1);
	}

	bool take(const T& t, size_t mult) {
		VERIFY(mult > 0);
		/*typename map<T, size_t>::iterator*/auto it = data_.find(t);
		if (it == data_.end()) {
			return false;
		} else {
			size_t have = it->second;
			if (have < mult) {
				data_.erase(it->first);
				return false;
			} else if (have == mult) {
				data_.erase(it->first);
				return true;
			} else {
				it->second -= mult;
				return true;
			}
		}
	}

	bool take(const T& t) {
		return take(t, 1);
	}

	size_t mult(const T& t) const {
		auto it = data_.find(t);
		if (it == data_.end()) {
			return 0;
		} else {
			return it->second;
		}
	}

	void clear() {
		data_.clear();
	}

	const_iterator begin() const {
		return data_.begin();
	}

	const_iterator end() const {
		return data_.end();
	}

};

template<class Graph>
class AssemblyPathCallback: public PathProcessor<Graph>::Callback {
public:
	typedef typename Graph::EdgeId EdgeId;

private:
	const Graph& g_;
	const ColorHandler<Graph>& coloring_;
	const edge_type assembly_color_;
	size_t edge_count_;

	std::vector<vector<EdgeId>> paths_;

	bool CheckPath(const vector<EdgeId>& path) const {
		if (path.size() > edge_count_)
			return false;
		for (auto it = path.begin(); it != path.end(); ++it) {
			if ((coloring_.Color(*it) & assembly_color_) == 0) {
				return false;
			}
		}
		return true;
	}

public:

	AssemblyPathCallback(const Graph& g, const ColorHandler<Graph>& coloring,
			edge_type assembly_color, size_t edge_count) :
			g_(g), coloring_(coloring), assembly_color_(assembly_color), edge_count_(edge_count) {
	}

	virtual void HandlePath(const vector<EdgeId>& path) {
		if (CheckPath(path)) {
			paths_.push_back(path);
		}
	}

	size_t size() const {
		return paths_.size();
	}

	std::vector<vector<EdgeId>> paths() const {
		return paths_;
	}
};

template<class Graph>
class SimpleInDelCorrector {
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	Graph& g_;
	ColorHandler<Graph>& coloring_;

	//become invalidated during process
//	const EdgesPositionHandler<Graph>& edge_pos_;
	vector<EdgeId> genome_path_;
	const edge_type genome_color_;
	const edge_type assembly_color_;

	bag<EdgeId> genome_edge_mult_;

	vector<EdgeId> FindAssemblyPath(VertexId start, VertexId end,
			size_t edge_count_bound, size_t min_length, size_t max_length) {
		AssemblyPathCallback<Graph> assembly_callback(g_, coloring_, assembly_color_,
				edge_count_bound);
		PathProcessor<Graph> path_finder(g_, min_length, max_length, start, end,
				assembly_callback);
		path_finder.Process();
		if (assembly_callback.size() == 1) {
			return assembly_callback.paths().front();
		}
		return {};
	}

	int TryFindGenomePath(size_t pos, VertexId end,
			size_t edge_count_bound) {
		for (size_t i = 0;
				i + pos < genome_path_.size() && i < edge_count_bound; ++i) {
			if (g_.EdgeEnd(genome_path_[pos + i]) == end) {
				return pos + i + 1;
			}
		}
		return -1;
	}

//	bag<edge_type> ColorLengths(const vector<EdgeId>& edges) {
//		bag<edge_type> answer;
//		for (size_t i = 0; i < edges.size(); ++i) {
//			answer.put(coloring_.Color(edges[i]), g_.length(edges[i]));
//		}
//		return answer;
//	}

	size_t VioletLengthOfGenomeUnique(const vector<EdgeId>& edges) {
		size_t answer = 0;
		for (size_t i = 0; i < edges.size(); ++i) {
			if (coloring_.Color(edges[i]) == edge_type::violet
					&& genome_edge_mult_.mult(edges[i]) == 1) {
				answer += g_.length(edges[i]);
			}
		}
		return answer;
	}

	//genome pos exclusive
//	size_t CumulativeGenomeLengthToPos(size_t pos) {
//		size_t answer = 0;
//		for (size_t i = 0; i < pos; ++i) {
//			answer += g_.length(genome_path_[i]);
//		}
//		return answer;
//	}

	bool CheckGenomePath(size_t genome_start, size_t genome_end) {
		return VioletLengthOfGenomeUnique(vector<EdgeId>(genome_path_.begin() + genome_start, genome_path_.begin() + genome_end)) < 25;
	}

	optional<pair<size_t, size_t>> FindGenomePath(VertexId start,
			VertexId end, size_t edge_count_bound) {
		for (size_t i = 0; i < genome_path_.size(); ++i) {
			if (g_.EdgeStart(genome_path_[i]) == start) {
				int path_end = TryFindGenomePath(i, end,
						edge_count_bound);
				if (path_end > 0 && CheckGenomePath(i, path_end))
					return make_optional(make_pair(size_t(i), size_t(path_end)));
			}
		}
		return boost::none;
	}

	void RemoveObsoleteEdges(const vector<EdgeId>& edges) {
		for (auto it = edges.begin(); it != edges.end(); ++it) {
			if (coloring_.Color(*it) == genome_color_ && genome_edge_mult_.mult(*it) == 0) {
				g_.DeleteEdge(*it);
				g_.CompressVertex(g_.EdgeStart(*it));
				if (!g_.RelatedVertices(g_.EdgeStart(*it), g_.EdgeEnd(*it))) {
					g_.CompressVertex(g_.EdgeEnd(*it));
				}
			}
		}
	}

	string GenomePathStr(size_t genome_start, size_t genome_end) const {
		return g_.str(vector<EdgeId>(genome_path_.begin() + genome_start,
				genome_path_.begin() + genome_end));
	}

	void CorrectGenomePath(size_t genome_start, size_t genome_end,
			const vector<EdgeId>& assembly_path) {
		DEBUG("Substituting genome path " << GenomePathStr(genome_start, genome_end) << " with assembly path " << g_.str(assembly_path))
		vector<EdgeId> genomic_edges;
		for (size_t i = genome_start; i < genome_end; ++i) {
			//side effects
			VERIFY(genome_edge_mult_.take(genome_path_[i]));
			genomic_edges.push_back(genome_path_[i]);
		}
		for (size_t i = 0; i < assembly_path.size(); ++i) {
			genome_edge_mult_.put(assembly_path[i]);
			coloring_.Paint(assembly_path[i], genome_color_);
		}
		genome_path_.insert(
				genome_path_.erase(genome_path_.begin() + genome_start,
						genome_path_.begin() + genome_end),
				assembly_path.begin(), assembly_path.end());
		RemoveObsoleteEdges(genomic_edges);
	}

//	pair<string, pair<size_t, size_t>> ContigIdAndPositions(EdgeId e) {
//		vector<EdgePosition> poss = edge_pos_.GetEdgePositions(e);
//		VERIFY(!poss.empty());
//		if (poss.size() > 1) {
//			WARN("Something strange with assembly positions");
//			return make_pair("", make_pair(0, 0));
//		}
//		EdgePosition pos = poss.front();
//		return make_pair(pos.contigId_, make_pair(pos.start(), pos.end()));
//	}

//	void WriteAltPath(EdgeId e, const vector<EdgeId>& genome_path) {
//		LengthIdGraphLabeler<Graph> basic_labeler(g_);
//		EdgePosGraphLabeler<Graph> pos_labeler(g_, edge_pos_);
//
//		CompositeLabeler<Graph> labeler(basic_labeler, pos_labeler);
//
//		string alt_path_folder = folder_ + ToString(g_.int_id(e)) + "/";
//		make_dir(alt_path_folder);
//		WriteComponentsAlongPath(g_, labeler, alt_path_folder + "path.dot", /*split_length*/
//		1000, /*vertex_number*/15, TrivialMappingPath(g_, genome_path),
//				*ConstructBorderColorer(g_, coloring_));
//	}

	//todo use contig constraints here!!!
	void AnalyzeGenomeEdge(EdgeId e) {
		DEBUG("Analysing shortcut genome edge " << g_.str(e));
		DEBUG("Multiplicity " << genome_edge_mult_.mult(e));
		if (genome_edge_mult_.mult(e) == 1) {
			vector<EdgeId> assembly_path = FindAssemblyPath(g_.EdgeStart(e), g_.EdgeEnd(e),
						100, 0, g_.length(e) + 1000);
			if (!assembly_path.empty()) {
				DEBUG("Assembly path " << g_.str(assembly_path));
				auto it = std::find(genome_path_.begin(), genome_path_.end(), e);
				VERIFY(it != genome_path_.end());
				size_t pos = it - genome_path_.begin();
				CorrectGenomePath(pos, pos + 1, assembly_path);
			} else {
				DEBUG("Couldn't find assembly path");
			}
		}
	}

	void AnalyzeAssemblyEdge(EdgeId e) {
		DEBUG("Analysing shortcut assembly edge " << g_.str(e));
		optional<pair<size_t, size_t>> genome_path = FindGenomePath(
				g_.EdgeStart(e), g_.EdgeEnd(e), /*edge count bound*//*100*/300);
		if (genome_path) {
			CorrectGenomePath(genome_path->first, genome_path->second,
					vector<EdgeId>{e});
		} else {
			DEBUG("Empty genome path");
		}
	}

	void FillGenomeEdgeMult() {
		for (auto it = genome_path_.begin(); it != genome_path_.end(); ++it) {
			genome_edge_mult_.put(*it);
		}
	}

public:
	SimpleInDelCorrector(Graph& g, ColorHandler<Graph>& coloring,
			const vector<EdgeId> genome_path,
			edge_type genome_color, edge_type assembly_color) :
			g_(g), coloring_(coloring), genome_path_(
					genome_path), genome_color_(genome_color), assembly_color_(assembly_color) {
	}

	void Analyze() {
		FillGenomeEdgeMult();
		for (auto it = g_.SmartEdgeBegin(); !it.IsEnd(); ++it) {
			if (coloring_.Color(*it) == genome_color_) {
				AnalyzeGenomeEdge(*it);
			}
			if (coloring_.Color(*it) == assembly_color_) {
				AnalyzeAssemblyEdge(*it);
			}
		}
	}

private:
	DECL_LOGGER("SimpleInDelCorrector")
	;
};

}
