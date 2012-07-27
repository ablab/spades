#pragma once

#include "utils.hpp";

namespace cap {

template<class Graph>
class AssemblyPathCallback: public PathProcessor<Graph>::Callback {
public:
	typedef typename Graph::EdgeId EdgeId;

private:
	const Graph& g_;
	size_t edge_count_;

	std::vector<vector<EdgeId>> paths_;

	bool CheckPath(const vector<EdgeId>& path) {

	}

public:

	AssemblyPathCallback(const Graph& g, size_t edge_count) :
			g_(g), edge_count_(edge_count) {
	}

	virtual void HandlePath(const vector<EdgeId>& path) {
		CheckPath(path);
		paths_.push_back(path);
	}

	size_t size() {
		return paths_.size();
	}

	std::vector<vector<EdgeId>> paths() {
		return paths_;
	}
};

template<class Graph>
class SimpleInDelCorrector {
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	const Graph& g_;
	const ColorHandler<Graph>& coloring_;

	//become invalidated during process
	const EdgesPositionHandler<Graph>& edge_pos_;
	const vector<EdgeId> genome_path_;
	const edge_type genome_color_;
	const edge_type assembly_color_;

	vector<EdgeId> FindAssemblyPath(VertexId start, VertexId end, size_t edge_count_bound
			, size_t min_length, size_t max_length) {
		PathProcessor<Graph> path_finder(g_, min_length, max_length, start, end);


	}

	vector<EdgeId> TryFindGenomePath(size_t pos, VertexId end,
			size_t edge_count_bound) {
		vector<EdgeId> answer;
		for (size_t i = 0;
				i + pos < genome_path_.size() && i < edge_count_bound; ++i) {
//			if ((coloring_.Color(genome_path_[pos + i]) & shortcut_color_) > 0) {
//				DEBUG("Came into edge of wrong color");
//				return vector<EdgeId>();
//			}
			answer.push_back(genome_path_[pos + i]);
			if (g_.EdgeEnd(genome_path_[pos + i]) == end) {
				return answer;
			}
		}
		DEBUG("Edge bound reached");
		return vector<EdgeId>();
	}

	map<edge_type, size_t> ColorLengths(const vector<EdgeId>& edges) {
		map<edge_type, size_t> answer;
		for (size_t i = 0; i < edges.size(); ++i) {
			answer[coloring_.Color(edges[i])] += g_.length(edges[i]);
		}
		return answer;
	}

	size_t VioletLengthOfGenomeUnique(const vector<EdgeId>& edges) {
		size_t answer = 0;
		for (size_t i = 0; i < edges.size(); ++i) {
			if (coloring_.Color(edges[i]) == edge_type::violet
					&& std::count(genome_path_.begin(), genome_path_.end(), edges[i]) == 1) {
				answer += g_.length(edges[i]);
			}
		}
		return answer;
	}

	//genome pos exclusive
	size_t CumulativeGenomeLengthToPos(size_t pos) {
		size_t answer = 0;
		for (size_t i = 0; i < pos; ++i) {
			answer += g_.length(genome_path_[i]);
		}
		return answer;
	}

	pair<vector<EdgeId>, pair<size_t, size_t>> FindGenomePath(VertexId start,
			VertexId end, size_t edge_count_bound) {
		for (size_t i = 0; i < genome_path_.size(); ++i) {
			if (g_.EdgeStart(genome_path_[i]) == start) {
				vector<EdgeId> path = TryFindGenomePath(i, end, edge_count_bound);
				if (!path.empty())
					return make_pair(
							path,
							make_pair(
									CumulativeGenomeLengthToPos(i),
									CumulativeGenomeLengthToPos(
											i + path.size())));
			}
		}
		return make_pair(vector<EdgeId>(), make_pair(0, 0));
	}

	pair<string, pair<size_t, size_t>> ContigIdAndPositions(EdgeId e) {
		vector<EdgePosition> poss = edge_pos_.GetEdgePositions(e);
		VERIFY(!poss.empty());
		if (poss.size() > 1) {
			WARN("Something strange with assembly positions");
			return make_pair("", make_pair(0, 0));
		}
		EdgePosition pos = poss.front();
		return make_pair(pos.contigId_, make_pair(pos.start(), pos.end()));
	}

	void WriteAltPath(EdgeId e, const vector<EdgeId>& genome_path) {
		LengthIdGraphLabeler<Graph> basic_labeler(g_);
		EdgePosGraphLabeler<Graph> pos_labeler(g_, edge_pos_);

		CompositeLabeler<Graph> labeler(basic_labeler, pos_labeler);

		string alt_path_folder = folder_ + ToString(g_.int_id(e)) + "/";
		make_dir(alt_path_folder);
		WriteComponentsAlongPath(g_, labeler, alt_path_folder + "path.dot", /*split_length*/
				1000, /*vertex_number*/15, TrivialMappingPath(g_, genome_path),
				*ConstructBorderColorer(g_, coloring_));
	}

	void Process(EdgeId e, const vector<EdgeId>& genome_path,
			size_t genome_start, size_t genome_end) {
		DEBUG("Processing edge and genome path");
		const size_t mem_lim = 2 << 26;
		Sequence edge_nucls = g_.EdgeNucls(e);
		Sequence path_nucls = MergeSequences(g_, genome_path);
		size_t edge_length = g_.length(e);
		size_t path_length = CummulativeLength(g_, genome_path);
		DEBUG(
				"Diff length " << abs((int) edge_length - (int) path_length)
						<< "; genome path length " << path_length
						<< "; edge length " << edge_length);
		pair<string, pair<size_t, size_t>> c_id_and_pos = ContigIdAndPositions(
				e);

		if (c_id_and_pos.first == "")
			return;

		WriteAltPath(e, genome_path);
		size_t unique_violet = VioletLengthOfGenomeUnique(genome_path);
		map<edge_type, size_t> color_cumm_lengths = ColorLengths(genome_path);
		if (color_cumm_lengths[edge_type::violet] * 10
				> color_cumm_lengths[edge_type::blue]) {
			WARN(
					"Very long path along violet edges: "
							<< color_cumm_lengths[edge_type::violet]
							<< " while blue path length: "
							<< color_cumm_lengths[edge_type::blue]);
			WARN("While processing edge: " << g_.str(e));
		}
		if (color_cumm_lengths[edge_type::violet] > 0)
			DEBUG("Violet edges in path");

		DEBUG("Total blue " << color_cumm_lengths[edge_type::blue]);
		DEBUG("Total violet " << color_cumm_lengths[edge_type::violet]);
		DEBUG("Total unique violet " << unique_violet);

		if (edge_length * path_length <= mem_lim) {
			size_t edit_dist = EditDistance(edge_nucls, path_nucls);
			DEBUG(
					"Edit distance " << edit_dist << ". That is "
							<< double(edit_dist) / max(edge_length, path_length));
			pair<size_t, size_t> local_sim = LocalSimilarity(edge_nucls,
					path_nucls);
			DEBUG(
					"Local sim " << local_sim.first << " interval length "
							<< local_sim.second << " relative "
							<< ((double) local_sim.first / local_sim.second));
//			assembly_length-genome_length relative_local_sim genome_path_length assembly_length genome_length
//			contig_id contig_start contig_end genome_start genome_end min max local_sim sim_interval edit_dist edit_dist/max tot_blue
//			tot_violet unique_violet edge_id
			cerr
					<< str(
							format(
									"%d %f %d %d %d %s %d %d %d %d %d %d %d %d %d %f %d %d %d %d")
									% ((int) edge_length - (int) path_length)
									% ((double) local_sim.first
											/ local_sim.second)
									% genome_path.size() % edge_length
									% path_length % c_id_and_pos.first
									% c_id_and_pos.second.first
									% c_id_and_pos.second.second % genome_start
									% genome_end % min(edge_length, path_length)
									% max(edge_length, path_length)
									% local_sim.first % local_sim.second
									% edit_dist
									% (double(edit_dist)
											/ max(edge_length, path_length))
									% color_cumm_lengths[edge_type::blue]
									% color_cumm_lengths[edge_type::violet]
									% unique_violet
									% g_.int_id(e)) << endl;
		} else {
			WARN("Edges were too long");
		}
	}

	void AnalyzeShortcutEdge(EdgeId e) {
		DEBUG("Analysing edge " << g_.str(e));
		pair<vector<EdgeId>, pair<size_t, size_t>> genome_path = FindGenomePath(
				g_.EdgeStart(e), g_.EdgeEnd(e), /*edge count bound*//*100*/300);
		if (!genome_path.first.empty()) {
			DEBUG(
					"Non empty genome path of edge count "
							<< genome_path.first.size());
			DEBUG("Path " << g_.str(genome_path.first));
			Process(e, genome_path.first, genome_path.second.first,
					genome_path.second.second);
		} else {
			DEBUG("Empty genome path");
		}
	}

public:
	SimpleInDelAnalyzer(const Graph& g, const ColorHandler<Graph>& coloring,
			const EdgesPositionHandler<Graph>& edge_pos,
			const vector<EdgeId> genome_path, edge_type shortcut_color,
			const string& folder) :
			g_(g), coloring_(coloring), edge_pos_(edge_pos), genome_path_(
					genome_path), shortcut_color_(shortcut_color), folder_(
					folder) {
	}

	void Analyze() {
		cerr
				<< "assembly_length-genome_length relative_local_sim genome_path_length assembly_length genome_length "
				<< "contig_id contig_start contig_end genome_start genome_end min max local_sim sim_interval edit_dist "
				<< "edit_dist/max tot_blue tot_violet unique_violet edge_id" << endl;
		for (auto it = g_.SmartEdgeBegin(); !it.IsEnd(); ++it) {
			if (coloring_.Color(*it) == shortcut_color_) {
				AnalyzeShortcutEdge(*it);
			}
		}
	}
private:
	DECL_LOGGER("SimpleInDelAnalyzer")
	;
};

}
