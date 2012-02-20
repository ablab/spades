#pragma once

#include "standard.hpp"
#include "utils.hpp"
#include "graph_pack.hpp"
#include "graph_construction.hpp"
#include "graph_simplification.hpp"
#include "simple_tools.hpp"
#include "omni/omni_utils.hpp"
#include "debruijn_stats.hpp"

namespace debruijn_graph {
typedef io::IReader<io::SingleRead> ContigStream;
typedef io::MultifileReader<io::SingleRead> CompositeContigStream;

enum edge_type {
	//don't change order!!!
	black = 0,
	red,
	blue,
	violet
};

template<class Graph, class Element>
class ElementColorHandler: public GraphActionHandler<Graph> {
	typedef GraphActionHandler<Graph> base;

	map<Element, edge_type> data_;
public:
	static string color_str(edge_type color) {
		static string colors[] = { "black", "red", "blue", "purple" };
		return colors[(int) color];
	}

	ElementColorHandler(const Graph& g) :
			base(g, "ElementColorHandler") {

	}

	void Paint(Element e, edge_type color) {
		data_[e] = (edge_type) ((int) data_[e] | (int) color);
	}

	edge_type Color(Element e) const {
		auto it = data_.find(e);
		if (it == data_.end())
			return edge_type::black;
		else
			return it->second;
	}

	string ColorStr(Element e) const {
		return color_str(Color(e));
	}

	/*virtual*/
	void HandleDelete(Element e) {
		data_.erase(e);
	}

//	virtual void HandleMerge(vector<EdgeId> old_edges, EdgeId new_edge) {
//		for(auto it = old_edges.begin(); it != old_edges.end(); ++it) {
//			PaintEdge(new_edge, Color(*it));
//		}
//	}
//
//	virtual void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) {
//		PaintEdge(new_edge, Color(edge1));
//		PaintEdge(new_edge, Color(edge2));
//	}
//
//	virtual void HandleSplit(EdgeId old_edge, EdgeId new_edge_1,
//			EdgeId new_edge2) {
//		PaintEdge(new_edge_1, Color(old_edge));
//		PaintEdge(new_edge2, Color(old_edge));
//	}

};

template<class Graph>
class ColorHandler: public GraphActionHandler<Graph> {
	typedef GraphActionHandler<Graph> base;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	ElementColorHandler<Graph, EdgeId> edge_color_;
	ElementColorHandler<Graph, VertexId> vertex_color_;
public:

	ColorHandler(const Graph& g) :
			base(g, "ColorHandler"), edge_color_(g), vertex_color_(g) {

	}

	void Paint(EdgeId e, edge_type color) {
		edge_color_.Paint(e, color);
	}

	void Paint(VertexId v, edge_type color) {
		vertex_color_.Paint(v, color);
	}

	edge_type Color(EdgeId e) const {
		return edge_color_.Color(e);
	}

	edge_type Color(VertexId v) const {
		return vertex_color_.Color(v);
	}

	map<EdgeId, string> EdgeColorMap() const {
		map<EdgeId, string> answer;
		for (auto it = this->g().SmartEdgeBegin(); !it.IsEnd(); ++it) {
			answer[*it] = edge_color_.ColorStr(*it);
		}
		return answer;
	}

	map<VertexId, string> VertexColorMap() const {
		map<VertexId, string> answer;
		for (auto it = this->g().begin(); it != this->g().end(); ++it) {
			answer[*it] = vertex_color_.ColorStr(*it);
		}
		return answer;
	}

};

template<class gp_t>
class CoveredRangesFinder {
	const gp_t& gp_;

	typedef typename gp_t::graph_t Graph;
	typedef typename Graph::EdgeId EdgeId;
	typedef map<EdgeId, vector<Range>> CoveredRanges;

	vector<Range> ProcessRange(Range new_range,
			const vector<Range>& curr_ranges) {
		vector<Range> answer;
		size_t i = 0;
		while (i < curr_ranges.size()
				&& curr_ranges[i].end_pos < new_range.start_pos) {
			answer.push_back(curr_ranges[i]);
			++i;
		}

		size_t merge_start =
				(i != curr_ranges.size()) ?
						std::min(curr_ranges[i].start_pos,
								new_range.start_pos) :
						new_range.start_pos;

		size_t merge_end = new_range.end_pos;
		while (i < curr_ranges.size()
				&& curr_ranges[i].start_pos <= new_range.end_pos) {
			if (curr_ranges[i].end_pos > merge_end)
				merge_end = curr_ranges[i].end_pos;
			++i;
		}
		answer.push_back(Range(merge_start, merge_end));
		while (i < curr_ranges.size()) {
			answer.push_back(curr_ranges[i]);
			++i;
		}
		return answer;
	}

	void ProcessPath(const MappingPath<EdgeId>& path, CoveredRanges& crs) {
		for (size_t i = 0; i < path.size(); ++i) {
			auto mapping = path[i];
			EdgeId edge = mapping.first;
			const vector<Range>& curr_ranges = crs[edge];
			Range mapping_range = mapping.second.mapped_range;
			VERIFY(gp_.g.length(edge) >= mapping_range.end_pos);
			crs[edge] = ProcessRange(mapping_range, curr_ranges);
			VERIFY(gp_.g.length(edge) >= crs[edge].back().end_pos);
		}
	}

public:

	CoveredRangesFinder(const gp_t& gp) :
			gp_(gp) {

	}

	void FindCoveredRanges(CoveredRanges& crs, ContigStream& stream) {
		io::SingleRead read;
		stream.reset();
		NewExtendedSequenceMapper<gp_t::k_value + 1, Graph> mapper(gp_.g,
				gp_.index, gp_.kmer_mapper);
		while (!stream.eof()) {
			stream >> read;
			ProcessPath(mapper.MapSequence(read.sequence()), crs);
		}
	}

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
class ComponentClassifier {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

public:
	enum component_type {
		error = 0,
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
		const string names[] = { "error", "single_red", "single_blue",
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
			g_(g), coloring_(coloring), bulge_length_(g_.k() * 5) {
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
		if (CheckTip(component)) {
			return component_type::tip;
		}
		if (CheckSimpleMisassembly(component))
			return component_type::simple_misassembly;
		if (CheckMonochrome(component))
			return component_type::monochrome;
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

	size_t to_draw_;
	ComponentClassifier<Graph> cc_;

public:
	ComponentTypeFilter(const Graph &g, size_t to_draw,
			const ColorHandler<Graph> &coloring) :
			base(g), to_draw_(to_draw), cc_(g, coloring) {
	}

	/*virtual*/
	bool Check(const vector<VertexId>& component) const {
		return cc_.GetComponentType(component) == to_draw_;
	}

private:
	DECL_LOGGER("ComponentTypeFilter")
	;
};

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
		make_dir("assembly_comparison");
		LongEdgesExclusiveSplitter<Graph> splitter(this->graph(), 1000000000);
		WriteComponents(this->graph(), splitter, *this,
				"breakpoint_graph", "assembly_comparison/breakpoint_graph.dot",
				coloring_.EdgeColorMap(), labeler);
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
	const Graph &graph_;
	const ColorHandler<Graph> &coloring_;
public:
	BPGraphStatCounter(const Graph &g, const ColorHandler<Graph> &coloring) :
			graph_(g), coloring_(coloring) {
	}

	void PrintComponents(size_t c_type,
			const GraphLabeler<Graph>& labeler) const {
		make_dir("assembly_comparison/");
		string type_dir = "assembly_comparison/"
				+ ComponentClassifier<Graph>::info_printer_pos_name(c_type)
				+ "/";
		make_dir(type_dir);
		string picture_dir = type_dir + "pictures/";
		make_dir(picture_dir);
		LongEdgesExclusiveSplitter<Graph> splitter(graph_, 1000000000);
		ComponentTypeFilter<Graph> stats(graph_, c_type, coloring_);
		WriteComponents(this->graph_, splitter, stats, "breakpoint_graph",
				picture_dir + "breakpoint_graph.dot", coloring_.EdgeColorMap(),
				labeler);
	}

	void PrintStats(const BreakPointGraphStatistics<Graph> &stats) const {
		make_dir("assembly_comparison/");
		for (size_t t = 0; t < component_type::size; t++) {
			string type_dir = "assembly_comparison/"
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

	void CountStats(const GraphLabeler<Graph>& labeler) const {
		make_dir("assembly_comparison/");
		BreakPointGraphStatistics<Graph> stats(graph_, coloring_);
		stats.CountStats();
		PrintStats(stats);
		PrintComponents(component_type::complex_misassembly, labeler);
		PrintComponents(component_type::monochrome, labeler);
		PrintComponents(component_type::tip, labeler);
		PrintComponents(component_type::simple_misassembly, labeler);
	}
};

template<class gp_t>
void FillPos(gp_t& gp, ContigStream& stream, string stream_prefix) {
	io::SingleRead read;
	while (!stream.eof()) {
		stream >> read;
		FillEdgesPos(gp, read.sequence(), stream_prefix + read.name());
	}
}

template<class gp_t>
class AssemblyComparer {
private:
	typedef typename gp_t::graph_t Graph;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef map<EdgeId, vector<Range>> CoveredRanges;
	typedef map<EdgeId, vector<size_t>> BreakPoints;

	void AddBreaks(set<size_t>& breaks, const vector<Range>& ranges) {
		for (auto it = ranges.begin(); it != ranges.end(); ++it) {
			breaks.insert(it->start_pos);
			breaks.insert(it->end_pos);
		}
	}

	vector<size_t> CombineCoveredRanges(const vector<Range>& ranges1,
			const vector<Range>& ranges2) {
		set<size_t> tmp_breaks;
		AddBreaks(tmp_breaks, ranges1);
		AddBreaks(tmp_breaks, ranges2);
		vector<size_t> breaks(tmp_breaks.begin(), tmp_breaks.end());
		//breaks contain 0 and edge_length here!
		VERIFY(breaks.size() >= 2);
		//cleaning breaks from 0 and edge_length
		vector<size_t> final_breaks;
		for (size_t i = 1; i < breaks.size() - 1; ++i) {
			final_breaks.push_back(breaks[i]);
		}
		return final_breaks;
	}

	void FindBreakPoints(const Graph& g, BreakPoints& bps, /*const */
			CoveredRanges& crs1, /*const */CoveredRanges& crs2) {
		for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
			EdgeId e = *it;
			bps[e] = CombineCoveredRanges(crs1[e], crs2[e]);
			VERIFY(bps[e].empty() || bps[e].back() < g.length(e));
		}
	}

	void SplitEdge(const vector<size_t>& breaks, EdgeId e, Graph& g) {
		vector<size_t> shifts(breaks.size());
		if (!breaks.empty()) {
			shifts[0] = breaks[0];
			for (size_t i = 1; i < breaks.size(); ++i) {
				shifts[i] = breaks[i] - breaks[i - 1];
			}
		}
		EdgeId curr_e = e;
		for (size_t i = 0; i < breaks.size(); ++i) {
			auto split_result = g.SplitEdge(curr_e, shifts[i]);
			curr_e = split_result.second;
		}
	}

	void SplitGraph(/*const */BreakPoints& bps, Graph& g) {
		set<EdgeId> initial_edges;
		for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
			initial_edges.insert(*it);
		}
		for (auto it = SmartSetIterator<Graph, EdgeId>(g, initial_edges.begin(),
				initial_edges.end()); !it.IsEnd(); ++it) {
			EdgeId e = *it;
			VERIFY(bps.find(e) != bps.end());
			VERIFY(bps[e].empty() || bps[e].back() < g.length(e));
			SplitEdge(bps[e], e, g);
		}
	}

	void DeleteVioletEdges(Graph& g, const ColorHandler<Graph>& coloring) {
		for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
			if (coloring.Color(*it) == edge_type::violet) {
				g.DeleteEdge(*it);
			}
		}
		Cleaner<Graph>(g).Clean();
	}

	void ColorPath(const Graph& g, const Path<EdgeId>& path, ColorHandler<Graph>& coloring,
			edge_type color) {
		for (size_t i = 0; i < path.size(); ++i) {
			coloring.Paint(path[i], color);
			coloring.Paint(g.EdgeStart(path[i]), color);
			coloring.Paint(g.EdgeEnd(path[i]), color);
		}
	}

	void ColorGraph(const gp_t& gp, ColorHandler<Graph>& coloring,
			ContigStream& stream, edge_type color) {
		io::SingleRead read;
		stream.reset();
		NewExtendedSequenceMapper<gp_t::k_value + 1, Graph> mapper(gp.g,
				gp.index, gp.kmer_mapper);
		while (!stream.eof()) {
			stream >> read;
			ColorPath(gp.g, mapper.MapSequence(read.sequence()).simple_path(),
					coloring, color);
		}
	}

public:

	void CompareAssemblies(ContigStream& stream1, ContigStream& stream2,
			const string& name1, const string& name2) {
		gp_t gp;
		CompositeContigStream stream(stream1, stream2);
		INFO("Constructing graph");
		ConstructGraph<gp_t::k_value, Graph>(gp.g, gp.index, stream);

//		debruijn_config::simplification::bulge_remover br_config;
//		br_config.max_bulge_length_coefficient = 4;
//		br_config.max_coverage = 1000.;
//		br_config.max_relative_coverage = 1.2;
//		br_config.max_delta = 3;
//		br_config.max_relative_delta = 0.1;
//		INFO("Removing bulges");
//		RemoveBulges(gp.g, br_config);

		INFO("Determining covered ranges");
		CoveredRangesFinder<gp_t> crs_finder(gp);
		CoveredRanges crs1;
		crs_finder.FindCoveredRanges(crs1, stream1);
		CoveredRanges crs2;
		crs_finder.FindCoveredRanges(crs2, stream2);
		BreakPoints bps;
		INFO("Determining breakpoints");
		FindBreakPoints(gp.g, bps, crs1, crs2);

		INFO("Splitting graph");
		SplitGraph(bps, gp.g);

		ColorHandler<Graph> coloring(gp.g);

		INFO("Coloring graph");
		ColorGraph(gp, coloring, stream1, edge_type::red);
		ColorGraph(gp, coloring, stream2, edge_type::blue);

		INFO("Filling contig positions");
		stream1.reset();
		FillPos(gp, stream1, name1);
		stream2.reset();
		FillPos(gp, stream2, name2);

		EdgePosGraphLabeler<Graph> pos_labeler(gp.g, gp.edge_pos);
		StrGraphLabeler<Graph> str_labeler(gp.g);
		CompositeLabeler<Graph> labeler(pos_labeler, str_labeler);

//		WriteToDotFile(gp.g, labeler, "oppa.dot", "oppa");

		LongEdgesInclusiveSplitter<Graph> splitter(gp.g, 10000);
		ComponentSizeFilter<Graph> filter(gp.g, 1000000000, 2);
		make_dir("assembly_comparison");
		make_dir("assembly_comparison/initial_pics/");
		WriteComponents(gp.g, splitter, filter, "breakpoint_graph",
				"assembly_comparison/initial_pics/breakpoint_graph.dot",
				coloring.EdgeColorMap(), MapColorer<Graph, VertexId>(coloring.VertexColorMap()), labeler);

		INFO("Removing unnecessary edges");
		DeleteVioletEdges(gp.g, coloring);

//		ReliableSplitter<Graph> splitter(gp.g, /*max_size*/100, /*edge_length_bound*/5000);
//		BreakPointsFilter<Graph> filter(gp.g, coloring, 3);
		INFO("Counting stats, outputting pictures");
		BPGraphStatCounter<Graph> counter(gp.g, coloring);
		counter.CountStats(labeler);
	}
private:
	DECL_LOGGER("AssemblyComparer")
	;
};

}
