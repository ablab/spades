#pragma once

namespace compare {

template<class Graph>
class Component {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	vector<size_t> edge_lengths_;
	vector<VertexId> component_;
public:
	Component(const Graph &g, const vector<VertexId> &component) :
			component_(component) {
		for (size_t i = 0; i < component.size(); i++)
			for (size_t j = 0; j < component.size(); j++) {
				vector<EdgeId> edges = g.GetEdgesBetween(component[i],
						component[j]);
				for (auto it = edges.begin(); it != edges.end(); ++it) {
					edge_lengths_.push_back(g.length(*it));
				}
			}
		std::sort(edge_lengths_.rbegin(), edge_lengths_.rend());
	}

	bool operator<(const Component<Graph> &that) const {
		size_t i = 0;
		while (i < this->edge_lengths_.size() && i < that.edge_lengths_.size()
				&& this->edge_lengths_[i] == that.edge_lengths_[i])
			i++;
		if (i == that.edge_lengths_.size())
			return false;
		if (i == this->edge_lengths_.size())
			return true;
		return this->edge_lengths_[i] < that.edge_lengths_[i];
	}

	bool operator==(const Component<Graph> &that) const {
		if (this->edge_lengths_.size() != that.edge_lengths_.size())
			return false;
		for (size_t i = 0; i < this->edge_lengths_.size(); i++)
			if (this->edge_lengths_[i] != that.edge_lengths_[i])
				return false;
		return true;
	}

	const vector<size_t> &edge_lengths() const {
		return edge_lengths_;
	}

};

template<class Stream, class Graph>
Stream &operator<<(Stream &stream, const Component<Graph> &component) {
	const vector<size_t> &lengths = component.edge_lengths();
	for (size_t i = 0; i < lengths.size(); i++) {
		stream << lengths[i] << " ";
	}
	return stream;
}

template<class Graph>
class ComponentClassifier {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

public:
	enum component_type {
		any = 0,
		error,
		single_red,
		single_blue,
		simple_bulge,
		tip,
		simple_misassembly,
		monochrome,
		complex_misassembly,
		size
	};

	static string info_printer_pos_name(size_t t) {
		const string names[] = { "all", "error", "single_red", "single_blue",
				"simple_bulge", "tip", "simple_misassembly", "monochrome",
				"complex_misassembly", "size" };
		return names[t];
	}

private:
	const Graph &g_;
	const ColorHandler<Graph> &coloring_;
	const size_t bulge_length_;

public:
	ComponentClassifier(const Graph &g, const ColorHandler<Graph> &coloring,
			size_t bulge_length) :
			g_(g), coloring_(coloring), bulge_length_(bulge_length) {
	}

	ComponentClassifier(const Graph &g, const ColorHandler<Graph> &coloring) :
			g_(g), coloring_(coloring), bulge_length_(g_.k() * 1000000) {
	}

	edge_type GetColour(EdgeId edge) const {
		return coloring_.Color(edge);
	}

	bool CheckSimpleMisassembly(const vector<VertexId> &component) const {
		if (component.size() == 4) {
			for (size_t i = 0; i < 4; i++)
				for (size_t j = i + 1; j < 4; j++) {
					vector<VertexId> sources;
					sources.push_back(component[i]);
					sources.push_back(component[j]);
					vector<VertexId> sinks;
					for (size_t k = 0; k < 4; k++) {
						if (k != i && k != j)
							sinks.push_back(component[k]);
					}
					if (CheckSimpleMisassembly(sources, sinks)) {
						return true;
					}
				}
		}
		return false;
	}

	bool CheckSimpleMisassembly(const vector<VertexId>& sources,
			const vector<VertexId>& sinks) const {
		if (sources.size() != 2 || sinks.size() != 2)
			return false;
		for (size_t i = 0; i < sources.size(); i++)
			for (size_t j = 0; j < sinks.size(); j++) {
				if (g_.GetEdgesBetween(sources[i], sinks[j]).size() != 1) {
					return false;
				}
			}
		for (size_t i = 0; i < 2; i++) {
			if (g_.GetEdgesBetween(sources[i], sources[1 - i]).size() != 0)
				return false;
			if (g_.GetEdgesBetween(sinks[i], sinks[1 - i]).size() != 0)
				return false;
		}
		for (size_t i = 0; i < 2; i++) {
			if (GetColour(g_.GetEdgesBetween(sources[i], sinks[0])[0])
					== GetColour(g_.GetEdgesBetween(sources[i], sinks[1])[0]))
				return false;
			if (GetColour(g_.GetEdgesBetween(sources[0], sinks[i])[0])
					== GetColour(g_.GetEdgesBetween(sources[1], sinks[i])[0]))
				return false;
		}
		return true;
	}

	bool CheckIsolated(edge_type colour,
			const vector<VertexId> &component) const {
		if (component.size() != 2)
			return false;
		vector<EdgeId> edges01 = g_.GetEdgesBetween(component[0], component[1]);
		vector<EdgeId> edges10 = g_.GetEdgesBetween(component[1], component[0]);
		vector<EdgeId> edges;
		edges.insert(edges.end(), edges01.begin(), edges01.end());
		edges.insert(edges.end(), edges10.begin(), edges10.end());
		if (edges.size() != 1) {
			return false;
		}
		return GetColour(edges[0]) == colour;
	}

	bool CheckBulge(const vector<VertexId> &component) const {
		if (component.size() != 2)
			return false;
		vector<EdgeId> edges01 = g_.GetEdgesBetween(component[0], component[1]);
		vector<EdgeId> edges10 = g_.GetEdgesBetween(component[1], component[0]);
		vector<EdgeId> edges;
		edges.insert(edges.end(), edges01.begin(), edges01.end());
		edges.insert(edges.end(), edges10.begin(), edges10.end());
		return (edges01.size() == 0 || edges10.size() == 0) && edges.size() == 2
				&& g_.length(edges[0]) < bulge_length_
				&& g_.length(edges[1]) < bulge_length_;
	}

	size_t EdgeNumber(const vector<VertexId> &component) const {
		size_t result = 0;
		for (size_t i = 0; i < component.size(); i++)
			for (size_t j = 0; j < component.size(); j++) {
				result += g_.GetEdgesBetween(component[i], component[j]).size();
			}
		return result;
	}

	bool Connected(VertexId v1, VertexId v2) const {
		return g_.GetEdgesBetween(v1, v2).size() > 0;
	}

	bool CheckTip(const vector<VertexId> &component) const {
		if (component.size() != 3)
			return false;
		if (EdgeNumber(component) != 2)
			return false;
		for (size_t i = 0; i < 3; i++) {
			if (CheckFork(component[i], component[(i + 1) % 3],
					component[(i + 2) % 3]))
				return true;
		}
		return false;
	}

	bool CheckFork(VertexId base, VertexId tip1, VertexId tip2) const {
		return (Connected(base, tip1) && Connected(base, tip2))
				|| (Connected(tip1, base) && Connected(tip2, base));
	}

	bool CheckMonochrome(const vector<VertexId> &component) const {
		set<edge_type> colours;
		for (size_t i = 0; i < component.size(); i++)
			for (size_t j = 0; j < component.size(); j++) {
				vector<EdgeId> edges = g_.GetEdgesBetween(component[i],
						component[j]);
				for (auto it = edges.begin(); it != edges.end(); ++it) {
					colours.insert(GetColour(*it));
				}
			}
		return colours.size() == 1;
	}

	component_type GetComponentType(const vector<VertexId> &component) const {
		if (component.size() < 2)
			return component_type::error;
		if (component.size() == 2) {
			if (CheckIsolated(edge_type::red, component))
				return component_type::single_red;
			if (CheckIsolated(edge_type::blue, component))
				return component_type::single_blue;
			if (CheckBulge(component)) {
				return component_type::simple_bulge;
			}
			return component_type::complex_misassembly;
		}
		if (CheckMonochrome(component))
			return component_type::monochrome;
		if (CheckTip(component)) {
			return component_type::tip;
		}
		if (CheckSimpleMisassembly(component))
			return component_type::simple_misassembly;
		return component_type::complex_misassembly;
	}
};

template<class Graph>
class ComponentTypeFilter: public GraphComponentFilter<Graph> {
private:
	typedef GraphComponentFilter<Graph> base;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef typename ComponentClassifier<Graph>::component_type component_type;

	component_type to_draw_;
	ComponentClassifier<Graph> cc_;

public:
	ComponentTypeFilter(const Graph &g, component_type to_draw,
			const ColorHandler<Graph> &coloring) :
			base(g), to_draw_(to_draw), cc_(g, coloring) {
	}

	/*virtual*/
	bool Check(const vector<VertexId>& component) const {
		return to_draw_ == component_type::any
				|| cc_.GetComponentType(component) == to_draw_;
	}

private:
	DECL_LOGGER("ComponentTypeFilter")
	;
};

template<class Graph>
class BreakPointsFilter: public GraphComponentFilter<Graph> {
	typedef GraphComponentFilter<Graph> base;
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	const ColorHandler<Graph> coloring_;
	size_t color_threshold_;
public:
	BreakPointsFilter(const Graph& graph, const ColorHandler<Graph>& coloring,
			size_t color_threshold) :
			base(graph), coloring_(coloring), color_threshold_(color_threshold) {

	}

	bool MultiColored(const GraphComponent<Graph>& component) const {
		set<edge_type> colors;
		for (auto it = component.e_begin(); it != component.e_end(); ++it) {
			colors.insert(coloring_.Color(*it));
		}
		return colors.size() >= color_threshold_;
	}

	/*virtual*/
	//todo change to set or GraphComponent and add useful protected methods
	bool Check(const vector<VertexId> &component_veritces) const {
		GraphComponent<Graph> component(this->graph(),
				component_veritces.begin(), component_veritces.end());
		return component.v_size() > 2 && MultiColored(component);
	}

};

template<class Graph>
class BreakPointGraphStatistics: public GraphComponentFilter<Graph> {
private:
	typedef GraphComponentFilter<Graph> base;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef typename ComponentClassifier<Graph>::component_type component_type;

	const ColorHandler<Graph> &coloring_;
	ComponentClassifier<Graph> cc_;
	mutable vector<vector<Component<Graph>>> components_;
	mutable vector<size_t> total_red_;
	mutable vector<size_t> total_blue_;
	bool ready_;

	void UpdateTotalEdgeLength(component_type t, const vector<VertexId>& component) const {
		for(size_t i = 0; i < component.size(); i++) {
			for(size_t j = 0; j < component.size(); j++) {
				vector<EdgeId> edges = this->graph().GetEdgesBetween(component[i], component[j]);
				for(auto it = edges.begin(); it != edges.end(); ++it) {
					if(coloring_.Color(*it) == edge_type::red)
					total_red_[t] += this->graph().length(*it);
					if(coloring_.Color(*it) == edge_type::blue)
					total_blue_[t] += this->graph().length(*it);
				}
			}
		}
	}

	void UpdateStats(const vector<VertexId>& component) const {
		component_type t = cc_.GetComponentType(component);
		Component<Graph> c(this->graph(), component);
		components_[t].push_back(c);
		UpdateTotalEdgeLength(t, component);
	}

public:
	BreakPointGraphStatistics(const Graph &g, const ColorHandler<Graph> &coloring) :
	base(g), coloring_(coloring), cc_(g, coloring), components_(component_type::size), total_red_(component_type::size), total_blue_(component_type::size), ready_(
			false) {
	}

	/*virtual*/bool Check(const vector<VertexId>& component) const {
		UpdateStats(component);
		return false;
	}

	void CountStats() {
		EmptyGraphLabeler<Graph> labeler;
		make_dir("assembly_compare");
		LongEdgesExclusiveSplitter<Graph> splitter(this->graph(), 1000000000);
		WriteComponents(this->graph(), splitter, *this,
				"assembly_compare/breakpoint_graph.dot",
				*ConstructColorer(coloring_), labeler);
		ready_ = true;
		for (size_t i = 0; i < component_type::size; ++i) {
			INFO("Number of components of type " << ComponentClassifier<Graph>::info_printer_pos_name(i) << " is " << GetComponentNumber(i));
			INFO("Total length of red edges in these components is " << total_red_[i]);
			INFO("Total length of blue edges in these components is " << total_blue_[i]);
		}
	}

	size_t GetComponentNumber(size_t t) const {
		return components_[t].size();
	}

	const vector<Component<Graph>> &GetComponents(size_t t) const {
		std::sort(components_[t].rbegin(), components_[t].rend());
		return components_[t];
	}

private:
	DECL_LOGGER("BreakPointGraphStatistics");
};

template<class Graph>
class BPGraphStatCounter {
private:
	typedef typename ComponentClassifier<Graph>::component_type component_type;
	typedef typename Graph::VertexId VertexId;
	const Graph &graph_;
	const ColorHandler<Graph> &coloring_;
	const string output_folder_;
public:
	BPGraphStatCounter(const Graph &g, const ColorHandler<Graph> &coloring,
			const string& output_folder) :
			graph_(g), coloring_(coloring), output_folder_(output_folder) {
	}

	void PrintComponents(component_type c_type,
			const GraphLabeler<Graph>& labeler,
			bool create_subdir = true) const {
		string filename;
		if (create_subdir) {
			make_dir(output_folder_);
			string type_dir = output_folder_
					+ ComponentClassifier<Graph>::info_printer_pos_name(c_type)
					+ "/";
			make_dir(type_dir);
			string picture_dir = type_dir + "pictures/";
			make_dir(picture_dir);
			filename = picture_dir + "breakpoint_graph.dot";
		} else {
			filename = output_folder_
					+ ComponentClassifier<Graph>::info_printer_pos_name(c_type)
					+ ".dot";
		}
		LongEdgesExclusiveSplitter<Graph> splitter(graph_, 1000000000);
		ComponentTypeFilter<Graph> stats(graph_, c_type, coloring_);

		WriteComponents(this->graph_, splitter, stats, filename,
				*ConstructColorer(coloring_), labeler);
	}

	void PrintStats(const BreakPointGraphStatistics<Graph> &stats) const {
		make_dir(output_folder_);
		for (size_t t = 0; t < component_type::size; t++) {
			string type_dir = output_folder_
					+ ComponentClassifier<Graph>::info_printer_pos_name(t)
					+ "/";
			make_dir(type_dir);
			ofstream stream;
			stream.open((type_dir + "components.txt").c_str());
			const vector<Component<Graph>> &components = stats.GetComponents(t);
			for (auto it = components.begin(); it != components.end(); ++it) {
				stream << *it << endl;
			}
			stream.close();
		}
	}

	void CountStats(const GraphLabeler<Graph>& labeler, bool detailed_output =
			true) const {
		make_dir(output_folder_);
		BreakPointGraphStatistics<Graph> stats(graph_, coloring_);
		stats.CountStats();
		if (detailed_output) {
			PrintStats(stats);
			PrintComponents(component_type::complex_misassembly, labeler);
			PrintComponents(component_type::monochrome, labeler);
			PrintComponents(component_type::tip, labeler);
			PrintComponents(component_type::simple_misassembly, labeler);
		}
		PrintComponents(component_type::any, labeler, detailed_output);
	}
};

template<class Graph>
class TrivialBreakpointFinder: public AbstractFilter<
		vector<typename Graph::VertexId>> {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	struct bp_comp {
		bp_comp(const Graph& g, const ColorHandler<Graph>& coloring) :
				g_(g), coloring_(coloring) {
		}

		bool operator()(VertexId v1, VertexId v2) {
			return MaxRedBlueIncLength(v1) > MaxRedBlueIncLength(v2);
		}

		size_t MaxRedBlueIncLength(VertexId v) {
			vector<EdgeId> edges;
			insert_all(edges, g_.IncomingEdges(v));
			insert_all(edges, g_.OutgoingEdges(v));
			return MaxRedBlueLength(edges);
		}

	private:
		const Graph& g_;
		const ColorHandler<Graph>& coloring_;

		size_t MaxRedBlueLength(const vector<EdgeId> edges) {
			size_t max_length = 0;
			for (auto it = edges.begin(); it != edges.end(); ++it) {
				if (coloring_.Color(*it) == edge_type::blue
						|| coloring_.Color(*it) == edge_type::red) {
					if (g_.length(*it) > max_length) {
						max_length = g_.length(*it);
					}
				}
			}
			VERIFY(max_length > 0);
			return max_length;
		}
	};

	const Graph& g_;
	const ColorHandler<Graph>& coloring_;
	const EdgesPositionHandler<Graph>& pos_;

	void ReportBreakpoint(VertexId v, const string& folder,
			const string& prefix) {
		TRACE("Vertex " << g_.str(v) << " identified as breakpoint");
		LengthIdGraphLabeler<Graph> basic_labeler(g_);
		EdgePosGraphLabeler<Graph> pos_labeler(g_, pos_);

		CompositeLabeler<Graph> labeler(basic_labeler, pos_labeler);
		VERIFY(g_.OutgoingEdgeCount(v) > 0);
		EdgeId e = g_.OutgoingEdges(v).front();
		WriteComponentsAroundEdge(g_, e,
				folder + prefix + ToString(g_.int_id(v)) + "_loc.dot",
				*ConstructBorderColorer(g_, coloring_), labeler);
	}

	bool CheckEdges(const vector<EdgeId>& edges) const {
		set<edge_type> colors;
		for (auto it = edges.begin(); it != edges.end(); ++it) {
			colors.insert(coloring_.Color(*it));
		}
		return edges.size() == 2 && colors.count(edge_type::blue) == 1
				&& NotTips(edges);
	}

	bool IsTip(VertexId v) const {
		return g_.IncomingEdgeCount(v) + g_.OutgoingEdgeCount(v) == 1;
	}

	bool IsTip(EdgeId e) const {
		return (IsTip(g_.EdgeStart(e)) || IsTip(g_.EdgeEnd(e)))
				&& g_.length(e) < 200;
	}

	bool NotTips(const vector<EdgeId>& edges) const {
		for (auto it = edges.begin(); it != edges.end(); ++it)
			if (IsTip(*it))
				return false;
		return true;
	}

public:
	TrivialBreakpointFinder(const Graph& g, const ColorHandler<Graph>& coloring,
			const EdgesPositionHandler<Graph>& pos) :
			g_(g), coloring_(coloring), pos_(pos) {
	}

	void FindBreakPoints(const string& folder) {
		vector<VertexId> breakpoints;
		for (auto it = g_.SmartEdgeBegin(); !it.IsEnd(); ++it) {
			if (coloring_.Color(*it) == edge_type::red) {
				if (CheckEdges(g_.OutgoingEdges(g_.EdgeStart(*it))))
					breakpoints.push_back(g_.EdgeStart(*it));
//					ReportBreakpoint(g_.EdgeStart(*it));
				if (CheckEdges(g_.IncomingEdges(g_.EdgeEnd(*it))))
					breakpoints.push_back(g_.EdgeEnd(*it));
//					ReportBreakpoint(g_.EdgeEnd(*it));
			}
		}
		bp_comp comp(g_, coloring_);
		sort(breakpoints.begin(), breakpoints.end(), comp);
		for (size_t i = 0; i < breakpoints.size(); ++i) {
			ReportBreakpoint(
					breakpoints[i],
					folder,
					ToString(i) + "_"
							+ ToString(comp.MaxRedBlueIncLength(breakpoints[i]))
							+ "_");
		}
	}

	virtual bool Check(const vector<typename Graph::VertexId> &vertices) const {
		GraphComponent<Graph> component(g_, vertices.begin(), vertices.end());
		for (auto it = component.e_begin(); it != component.e_end(); ++it) {
			if (coloring_.Color(*it) == edge_type::red) {
				if (CheckEdges(g_.OutgoingEdges(g_.EdgeStart(*it)))
						|| CheckEdges(g_.IncomingEdges(g_.EdgeEnd(*it))))
					return true;
			}
		}
		return false;
	}

private:
	DECL_LOGGER("TrivialBreakpointFinder");
};

template<class Graph>
class SimpleInDelAnalyzer {
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	const Graph& g_;
	const ColorHandler<Graph>& coloring_;
	const vector<EdgeId> genome_path_;
	const edge_type shortcut_color_;

	vector<EdgeId> TryFindPath(size_t pos, VertexId end, size_t edge_count_bound) {
		vector<EdgeId> answer;
		for(size_t i = 0; i + pos < genome_path_.size() && i < edge_count_bound; ++i) {
			if ((coloring_.Color(genome_path_[pos + i]) & shortcut_color_) > 0)
				return vector<EdgeId>();
			answer.push_back(genome_path_[pos + i]);
			if (g_.EdgeEnd(genome_path_[pos + i]) == end) {
				return answer;
			}
		}
		return vector<EdgeId>();
	}

	vector<EdgeId> FindGenomePath(VertexId start, VertexId end, size_t edge_count_bound) {
		for (size_t i = 0; i < genome_path_.size(); ++i) {
			if (g_.EdgeStart(genome_path_[i]) == start) {
				vector<EdgeId> path = TryFindPath(i, end, edge_count_bound);
				if (!path.empty())
					return path;
			}
		}
		return vector<EdgeId>();
	}

	void Process(EdgeId e, const vector<EdgeId>& genome_path) {
		DEBUG("Processing edge and genome path");
		const size_t mem_lim = 2 << 26;
		Sequence edge_nucls = g_.EdgeNucls(e);
		Sequence path_nucls = MergeSequences(g_, genome_path);
		size_t edge_length = g_.length(e);
		size_t path_length = CummulativeLength(g_, genome_path);
		DEBUG("Diff length " << abs((int)edge_length - (int)path_length)
				<< "; genome path length " << path_length << "; edge length " << edge_length);
		if (edge_length * path_length <= mem_lim) {
			size_t edit_dist = EditDistance(edge_nucls, path_nucls);
			DEBUG("Edit distance " << edit_dist
					<< ". That is " << double(edit_dist) / max(edge_length, path_length));
			pair<size_t, size_t> local_sim = LocalSimilarity(edge_nucls, path_nucls);
			DEBUG("Local sim " << local_sim.first << " interval length " << local_sim.second << " relative " << ((double)local_sim.first/local_sim.second));
//			assembly_length-genome_length relative_local_sim genome_path_length assembly_length genome_length min max local_sim sim_interval edit_dist edit_dist/max
			cerr << str(format("%d %f %d %d %d %d %d %d %d %d %f")
					% ((int)edge_length - (int)path_length)
					% ((double)local_sim.first/local_sim.second)
					% genome_path.size()
					% edge_length
					% path_length
					% min(edge_length, path_length)
					% max(edge_length, path_length)
					% local_sim.first
					% local_sim.second
					% edit_dist
					% (double(edit_dist) / max(edge_length, path_length))) << endl;
		} else {
			WARN("Edges were too long");
		}
	}

	void AnalyzeShortcutEdge(EdgeId e) {
		DEBUG("Analysing edge " << g_.str(e));
		vector<EdgeId> genome_path = FindGenomePath(g_.EdgeStart(e), g_.EdgeEnd(e), /*edge count bound*/ 100);
		if (!genome_path.empty()) {
			DEBUG("Non empty genome path of edge count " << genome_path.size());
			DEBUG("Path " << g_.str(genome_path));
			Process(e, genome_path);
		} else {
			DEBUG("Empty genome path");
		}
	}

public:
	SimpleInDelAnalyzer(const Graph& g, const ColorHandler<Graph>& coloring,
			const vector<EdgeId> genome_path, edge_type shortcut_color)
	: g_(g), coloring_(coloring), genome_path_(genome_path), shortcut_color_(shortcut_color) {
	}

	void Analyze() {
		for (auto it = g_.SmartEdgeBegin(); !it.IsEnd(); ++it) {
			if (coloring_.Color(*it) == shortcut_color_) {
				AnalyzeShortcutEdge(*it);
			}
		}
	}
private:
	DECL_LOGGER("SimpleInDelAnalyzer");
};

template<class gp_t>
class SimpleRearrangementDetector {
private:
	typedef typename gp_t::graph_t Graph;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	const gp_t& gp_;
	const ColorHandler<Graph>& coloring_;
	const string ref_prefix_;
	const string folder_;
	mutable size_t count_;

	void ReportPossibleRearrangementConnection(EdgeId e, int start_ref_pos,
			int end_ref_pos, const string& folder) const {
		INFO(
				"Edge " << gp_.g.str(e)
						<< " identified as rearrangement connection");
		LengthIdGraphLabeler<Graph> basic_labeler(gp_.g);
		EdgePosGraphLabeler<Graph> pos_labeler(gp_.g, gp_.edge_pos);

		CompositeLabeler<Graph> labeler(basic_labeler, pos_labeler);

		INFO(
				count_ << " example start_ref_pos: " << start_ref_pos
						<< " end_ref_pos: " << end_ref_pos);
		string filename = str(
				boost::format("%s%d_%d_%d_%d.dot") % folder % count_
						% gp_.g.int_id(e) % start_ref_pos % end_ref_pos);
		WriteComponentsAroundEdge(gp_.g, e, filename,
				*ConstructBorderColorer(gp_.g, coloring_), labeler);
		count_++;
	}

	bool ContainsBlueEdge(const vector<EdgeId>& edges) const {
		for (size_t i = 0; i < edges.size(); ++i) {
			if (coloring_.Color(edges[i]) == edge_type::blue)
				return true;
		}
		return false;
	}

	EdgeId GetBlueEdge(const vector<EdgeId>& edges) const {
		for (size_t i = 0; i < edges.size(); ++i) {
			if (coloring_.Color(edges[i]) == edge_type::blue)
				return edges[i];
		}
		VERIFY(false);
		return EdgeId(NULL);
	}

	int GetRefPosition(EdgeId e, bool start_position) const {
		EdgePosition pos =
				RefPositions(gp_.edge_pos.GetEdgePositions(e)).front();
		int coeff = boost::ends_with(pos.contigId_, "_RC") ? -1 : 1;
		Range range = pos.m_range_.initial_range;
		return coeff * (start_position ? range.start_pos : range.end_pos);
	}

	bool IsSingleRefPosition(EdgeId e) const {
		return RefPositions(gp_.edge_pos.GetEdgePositions(e)).size() == 1;
	}

	vector<EdgePosition> RefPositions(const vector<EdgePosition>& poss) const {
		vector < EdgePosition > answer;
		for (auto it = poss.begin(); it != poss.end(); ++it) {
			if (boost::starts_with(it->contigId_, ref_prefix_)) {
				answer.push_back(*it);
			}
		}
		return answer;
	}

public:
	SimpleRearrangementDetector(const gp_t& gp,
			const ColorHandler<Graph>& coloring, const string& ref_prefix,
			const string& folder) :
			gp_(gp), coloring_(coloring), ref_prefix_(ref_prefix), folder_(
					folder), count_(0) {
	}

	void Detect() const {
		for (auto it = gp_.g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
			if (coloring_.Color(*it) == edge_type::red) {
				INFO("Processing red edge " << gp_.g.str(*it));
				if (gp_.g.OutgoingEdgeCount(gp_.g.EdgeStart(*it)) == 2
						&& ContainsBlueEdge(
								gp_.g.OutgoingEdges(gp_.g.EdgeStart(*it)))) {
					EdgeId first_edge = GetBlueEdge(
							gp_.g.OutgoingEdges(gp_.g.EdgeStart(*it)));
					if (gp_.g.IncomingEdgeCount(gp_.g.EdgeEnd(*it)) == 2
							&& ContainsBlueEdge(
									gp_.g.IncomingEdges(gp_.g.EdgeEnd(*it)))) {
						EdgeId second_edge = GetBlueEdge(
								gp_.g.IncomingEdges(gp_.g.EdgeEnd(*it)));
						if (first_edge != second_edge) {
							INFO("Edges passed topology checks");
							if (IsSingleRefPosition(first_edge)
									&& IsSingleRefPosition(second_edge)) {
								int start_ref_pos = GetRefPosition(first_edge,
										true);
								int end_ref_pos = GetRefPosition(second_edge,
										false);
								INFO("Edges had multiplicity one in reference");
								ReportPossibleRearrangementConnection(*it,
										start_ref_pos, end_ref_pos, folder_);
							} else {
								INFO("Ooops");
								INFO(
										"Edges had multiplicity more than one in reference");
							}
						}
					}
				}
			}
		}
	}

	DECL_LOGGER("SimpleRearrangementDetector");
};

template<class Graph>
class GraphEdgeEnumerator {
	const Graph& g_;
	typedef typename Graph::EdgeId EdgeId;
protected:
	GraphEdgeEnumerator(const Graph& g) :
			g_(g) {
	}

	const Graph& g() {
		return g_;
	}
public:
	virtual ~GraphEdgeEnumerator() {
	}
	virtual map<EdgeId, string> Enumerate() const;
};

template<class Graph>
class ThreadedGenomeEnumerator: public GraphEdgeEnumerator<Graph> {
	typedef GraphEdgeEnumerator<Graph> base;
	typedef typename Graph::EdgeId EdgeId;
	vector<EdgeId> genome_path_;
public:
	ThreadedGenomeEnumerator(const Graph& g, const vector<EdgeId>& genome_path) :
			base(g), genome_path_(genome_path) {
	}

	/*virtual */
	map<EdgeId, string> Enumerate() const {
		map < EdgeId, string > answer;
		//numerating genome path
		int curr = 0;
		for (auto it = genome_path_.begin(); it != genome_path_.end(); ++it) {
			if (answer.find(*it) == answer.end()) {
				curr++;
				answer[*it] = ToString(curr);
				answer[this->g().conjugate(*it)] = ToString(-curr);
			}
		}
		curr = 1000000;
		for (auto it = this->g().SmartEdgeBegin(); !it.IsEnd(); ++it) {
			if (answer.find(*it) == answer.end()) {
				curr++;
				answer[*it] = ToString(curr);
				answer[this->g().conjugate(*it)] = ToString(-curr);
			}
		}
		return answer;
	}
};

}
