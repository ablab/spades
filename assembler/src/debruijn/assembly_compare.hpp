#pragma once

#include "standard.hpp"
#include "utils.hpp"
#include "graph_pack.hpp"
#include "graph_construction.hpp"
#include "graph_simplification.hpp"
#include "simple_tools.hpp"
#include "omni/omni_utils.hpp"
#include "debruijn_stats.hpp"
#include "io/delegating_reader_wrapper.hpp"
#include "io/splitting_wrapper.hpp"

namespace debruijn_graph {
using namespace omnigraph;

template<class Graph>
Sequence MergeSequences(const Graph& g,
		const vector<typename Graph::EdgeId>& edges) {
	vector<const Sequence*> path_sequences;
	for (auto it = edges.begin(); it != edges.end(); ++it) {
		path_sequences.push_back(&g.EdgeNucls(*it));
	}
	return MergeOverlappingSequences(path_sequences, g.k());
}

template<class gpt>
class TipsProjector {
	typedef typename gpt::graph_t Graph;
	typedef typename Graph::EdgeId EdgeId;

	gpt& gp_;
	const UniquePathFinder<Graph> unique_path_finder_;

	optional<EdgeId> UniqueAlternativeEdge(EdgeId tip, bool outgoing_tip) {
		vector<EdgeId> edges;
		if (outgoing_tip) {
			edges = gp_.g.OutgoingEdges(gp_.g.EdgeStart(tip));
		} else {
			edges = gp_.g.IncomingEdges(gp_.g.EdgeEnd(tip));
		}
		set<EdgeId> edges_set(edges.begin(), edges.end());
		edges_set.erase(tip);
		if (edges_set.size() == 1)
			return optional<EdgeId>(*edges_set.begin());
		else
			return boost::none;
	}

	vector<EdgeId> UniqueAlternativePath(EdgeId tip, bool outgoing_tip) {
		optional<EdgeId> alt_edge = UniqueAlternativeEdge(tip, outgoing_tip);
		if (alt_edge) {
			if (outgoing_tip) {
				return unique_path_finder_.UniquePathForward(*alt_edge);
			} else {
				return unique_path_finder_.UniquePathBackward(*alt_edge);
			}
		}
		return vector<EdgeId>();
	}

	void AlignAndProject(const Sequence& tip_seq,
			const Sequence& alt_seq, bool outgoing_tip) {
		//todo refactor
		Sequence aligned_tip = tip_seq;
		Sequence aligned_alt = alt_seq;
		if (outgoing_tip) {
			if (tip_seq.size() >= alt_seq.size()) {
				aligned_tip = tip_seq.Subseq(0, alt_seq.size());
			} else {
				aligned_alt = alt_seq.Subseq(0, tip_seq.size());
			}
		} else {
			if (tip_seq.size() >= alt_seq.size()) {
				aligned_tip = tip_seq.Subseq(tip_seq.size() - alt_seq.size());
			} else {
				aligned_alt = alt_seq.Subseq(alt_seq.size() - tip_seq.size());
			}
		}
		gp_.kmer_mapper.RemapKmers(aligned_tip, aligned_alt);
	}

public:
	TipsProjector(gpt& gp) :
			gp_(gp), unique_path_finder_(gp.g) {

	}

	void ProjectTip(EdgeId tip) {
		bool outgoing_tip = gp_.g.IsDeadEnd(gp_.g.EdgeEnd(tip));
		Sequence tip_seq = gp_.g.EdgeNucls(tip);
		Sequence alt_seq = MergeSequences(gp_.g,
				UniqueAlternativePath(tip, outgoing_tip));
		if (tip_seq.size() > alt_seq.size()) {
			WARN("Can't fully project tip " << gp_.g.int_id(tip));
		} else {
			AlignAndProject(tip_seq, alt_seq, outgoing_tip);
		}
	}
private:
	DECL_LOGGER("TipsProjector")
	;
};


class ModifyingWrapper : public io::DelegatingReaderWrapper<io::SingleRead> {
public:
	ModifyingWrapper(io::IReader<io::SingleRead>& reader) : io::DelegatingReaderWrapper<io::SingleRead>(reader) {
	}

	virtual Sequence Modify(Sequence read) = 0;

	ModifyingWrapper& operator>>(io::SingleRead& read) {
		this->reader() >> read;
		Sequence s = read.sequence();
		read.SetSequence(Modify(s).str().c_str());
		return *this;
	}
};

template<size_t k, class Graph>
class ContigRefiner: public io::DelegatingReaderWrapper<io::SingleRead> {
	typedef io::DelegatingReaderWrapper<io::SingleRead> base;

	const Graph& graph_;
	NewExtendedSequenceMapper<k + 1, Graph> mapper_;

	Path<EdgeId> FixPath(const Path<EdgeId>& path) {
		for (size_t i = 1; i < path.size(); ++i) {
			VERIFY(graph_.EdgeStart(path[i]) == graph_.EdgeEnd(path[i - 1]));
		}
		return path;
	}

	Sequence Refine(const Sequence& s) {
		if(s < !s)
			return !Refine(!s);
		Path<EdgeId> path = FixPath(mapper_.MapSequence(s).simple_path());
		vector<EdgeId> edges = path.sequence();
		Sequence path_sequence = MergeSequences(graph_, edges);
		size_t start = path.start_pos();
		size_t end = path_sequence.size() - graph_.k()
				- graph_.length(edges.back()) + path.end_pos();
		return path_sequence.Subseq(start, end);
	}

public:
	ContigRefiner(io::IReader<io::SingleRead>& reader, const Graph &graph, const NewExtendedSequenceMapper<k + 1, Graph> &mapper) :
			base(reader), graph_(graph), mapper_(mapper) {
	}

	/* virtual */
	ContigRefiner& operator>>(io::SingleRead& read) {
		this->reader() >> read;
		Sequence s = Refine(read.sequence());
		read.SetSequence(s.str().c_str());
		string quality(s.size(), (char)33);
		read.SetQuality(quality.c_str());
		return *this;
	}

};

template<class gp_t>
class RefineGPConstructor {
private:
	typedef typename gp_t::graph_t Graph;
	typedef ContigRefiner<gp_t::k_value, Graph> Refiner;
	gp_t gp_;
public:
	RefineGPConstructor(io::IReader<io::SingleRead>& reader1, io::IReader<io::SingleRead>& reader2) {
		io::MultifileReader<io::SingleRead> composite_reader(reader1, reader2);
		ConstructGraph<gp_t::k_value>(gp_.g, gp_.index, composite_reader);

		//todo configure!!!
		debruijn_config::simplification::bulge_remover br_config;
		br_config.max_bulge_length_coefficient = 10;
		br_config.max_coverage = 1000.;
		br_config.max_relative_coverage = 1.2;
		br_config.max_delta = 5;
		br_config.max_relative_delta = 0.1;
		INFO("Removing bulges");
		RemoveBulges(gp_.g, br_config);

		TipsProjector<gp_t> tip_projector(gp_);
		boost::function<void(EdgeId)> projecting_callback =
				boost::bind(&TipsProjector<gp_t>::ProjectTip, &tip_projector, _1);
		debruijn_config::simplification::tip_clipper tc_config;
		tc_config.max_coverage = 1000.;
		tc_config.max_relative_coverage = 1.1;
		tc_config.max_tip_length_coefficient = 2.;

		ClipTips(gp_.g, tc_config, /*read_length*/100, projecting_callback);

		reader1.reset();
		reader2.reset();
	}

	const gp_t &GetGP() {
		return gp_;
	}

	NewExtendedSequenceMapper<gp_t::k_value + 1, Graph> ConstructMapper() {
		return NewExtendedSequenceMapper<gp_t::k_value + 1, Graph>(gp_.g, gp_.index, gp_.kmer_mapper);
	}
};

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
	typedef typename Graph::VertexId VertexId;
	const Graph &graph_;
	const ColorHandler<Graph> &coloring_;
public:
	BPGraphStatCounter(const Graph &g, const ColorHandler<Graph> &coloring) :
			graph_(g), coloring_(coloring) {
	}

	void PrintComponents(size_t c_type, const GraphLabeler<Graph>& labeler,
			const AbstractColorer<Graph, VertexId> &v_colorer) const {
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
				v_colorer, labeler);
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

	void CountStats(const GraphLabeler<Graph>& labeler,
			const AbstractColorer<Graph, VertexId> &v_colorer) const {
		make_dir("assembly_comparison/");
		BreakPointGraphStatistics<Graph> stats(graph_, coloring_);
		stats.CountStats();
		PrintStats(stats);
		PrintComponents(component_type::complex_misassembly, labeler,
				v_colorer);
		PrintComponents(component_type::monochrome, labeler, v_colorer);
		PrintComponents(component_type::tip, labeler, v_colorer);
		PrintComponents(component_type::simple_misassembly, labeler, v_colorer);
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

class RCSplittingStream : public io::DelegatingReaderWrapper<io::SingleRead> {
private:
	io::SplittingWrapper filtered_reader_;
	io::RCReaderWrapper<io::SingleRead> stream_;
public:
	RCSplittingStream(ContigStream &base_stream) : filtered_reader_(base_stream), stream_(filtered_reader_){
		Init(stream_);
	}
};

template<class gp_t>
class UntangledGraphConstructor {
private:
	typedef typename gp_t::graph_t Graph;
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

public:
	Graph graph_;
	IdTrackHandler<Graph> int_ids_;
	ColorHandler<Graph> untangled_coloring_;
	map<Sequence, vector<EdgeId>> red_paths_;
	map<Sequence, vector<EdgeId>> blue_paths_;
	map<EdgeId, EdgeId> purple_map_;
	map<VertexId, VertexId> vertex_map_;

private:
	VertexId GetStartVertex(const Graph &graph, const Path<EdgeId> &path, size_t i) {
		if(i != 0 || path.start_pos() == 0)
			return vertex_map_[graph.EdgeStart(path[i])];
		else
			return graph_.AddVertex();
	}

	VertexId GetEndVertex(const Graph &graph, const Path<EdgeId> &path, size_t i) {
		if(i != path.size() - 1 || path.end_pos() == graph.length(path[i]))
			return vertex_map_[graph.EdgeEnd(path[i])];
		else
			return graph_.AddVertex();
	}

	void Untangle(const gp_t& gp, const ColorHandler<Graph> &coloring, ContigStream& stream, edge_type color) {
		io::SingleRead read;
		stream.reset();
		while (!stream.eof()) {
			stream >> read;
			Untangle(gp, coloring, read.sequence(), color);
		}
	}

	void Untangle(const gp_t& gp, const ColorHandler<Graph> &coloring, Sequence contig, edge_type color) {
		NewExtendedSequenceMapper<gp_t::k_value + 1, Graph> mapper(gp.g, gp.index, gp.kmer_mapper);
		Path<EdgeId> path = mapper.MapSequence(contig).simple_path();
		vector<EdgeId> new_path;
		for(size_t i = 0; i < path.size(); i++) {
			EdgeId next;
			if(coloring.Color(path[i]) != edge_type::violet) {
				size_t j = i;
				vector<EdgeId> to_glue;
				while(j < path.size() && coloring.Color(path[j]) != edge_type::violet) {
					to_glue.push_back(path[j]);
					j++;
				}
				Sequence new_edge_sequence = MergeSequences(gp.g, to_glue);
				next = graph_.AddEdge(GetStartVertex(gp.g, path, i), GetEndVertex(gp.g, path, j - 1), new_edge_sequence);
				i = j - 1;
			} else {
				next = purple_map_[path[i]];
			}
			new_path.push_back(next);
			untangled_coloring_.Paint(next, color);
			untangled_coloring_.Paint(graph_.EdgeStart(next), color);
			untangled_coloring_.Paint(graph_.EdgeEnd(next), color);
		}
		if(color == edge_type::red)
			red_paths_[contig] = new_path;
		else
			blue_paths_[contig] = new_path;
	}

public:
	UntangledGraphConstructor(const gp_t &gp,
			const ColorHandler<Graph> &coloring,
			io::IReader<io::SingleRead> &stream1,
			io::IReader<io::SingleRead> &stream2) :
			graph_(gp_t::k_value), int_ids_(graph_),
			untangled_coloring_(graph_) {
		const Graph & old_graph = gp.g;
		for (auto it = old_graph.begin(); it != old_graph.end(); ++it) {
			vertex_map_[*it] = graph_.AddVertex();
		}
		for (auto it = old_graph.SmartEdgeBegin(); !it.IsEnd(); ++it) {
			if (coloring.Color(*it) == edge_type::violet) {
				EdgeId new_edge = graph_.AddEdge(
						vertex_map_[old_graph.EdgeStart(*it)],
						vertex_map_[old_graph.EdgeEnd(*it)],
						old_graph.EdgeNucls(*it));
				purple_map_[*it] = new_edge;
				untangled_coloring_.Paint(new_edge, edge_type::violet);
				untangled_coloring_.Paint(graph_.EdgeStart(new_edge), edge_type::violet);
				untangled_coloring_.Paint(graph_.EdgeEnd(new_edge), edge_type::violet);
			}
		}
		Untangle(gp, coloring, stream1, edge_type::red);
		Untangle(gp, coloring, stream2, edge_type::blue);
	}

};

template<class gp_t>
class AssemblyComparer {
private:
	typedef typename gp_t::graph_t Graph;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef map<EdgeId, vector<Range>> CoveredRanges;
	typedef map<EdgeId, vector<size_t>> BreakPoints;

	io::IReader<io::SingleRead> &stream1_;
	io::IReader<io::SingleRead> &stream2_;

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

	void SplitGraph(gp_t& gp) {
		INFO("Determining covered ranges");
		CoveredRangesFinder<gp_t> crs_finder(gp);
		CoveredRanges crs1;
		crs_finder.FindCoveredRanges(crs1, stream1_);
		CoveredRanges crs2;
		crs_finder.FindCoveredRanges(crs2, stream2_);
		BreakPoints bps;
		INFO("Determining breakpoints");
		FindBreakPoints(gp.g, bps, crs1, crs2);

		INFO("Splitting graph");
		SplitGraph(bps, gp.g);
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

	void ColorPath(const Graph& g, const Path<EdgeId>& path,
			ColorHandler<Graph>& coloring, edge_type color) {
		for (size_t i = 0; i < path.size(); ++i) {
			coloring.Paint(path[i], color);
			coloring.Paint(g.EdgeStart(path[i]), color);
			coloring.Paint(g.EdgeEnd(path[i]), color);
		}
	}

	void ColorGraph(const gp_t& gp, ColorHandler<Graph>& coloring) {
		INFO("Coloring graph");
		ColorGraph(gp, coloring, stream1_, edge_type::red);
		ColorGraph(gp, coloring, stream2_, edge_type::blue);
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

	void SimplifyGraph(Graph & graph) {
		debruijn_config::simplification::bulge_remover br_config;
		br_config.max_bulge_length_coefficient = 50;
		br_config.max_coverage = 1000.;
		br_config.max_relative_coverage = 1.2;
		br_config.max_delta = 20;
		br_config.max_relative_delta = 0.1;
		INFO("Removing bulges");
		RemoveBulges(graph, br_config);

//		debruijn_config::simplification::tip_clipper tc;
//		tc.max_coverage = 1000;
//		tc.max_relative_coverage = 1000;
//		tc.max_tip_length_coefficient = 6;
//		ClipTips(gp.g, tc, 10 * gp.g.k());
	}

	void SaveGraph(string path, gp_t &gp) {
		INFO("Saving graph");
		PrintGraphPack("assembly_comparison/saves/graph", gp);
	}

public:

	AssemblyComparer(io::IReader<io::SingleRead> &raw_reader_1,
			io::IReader<io::SingleRead> &raw_reader_2) : stream1_(raw_reader_1), stream2_(raw_reader_2) {
	}

	void CompareAssemblies(const string& name1, const string& name2) {
		stream1_.reset();
		stream2_.reset();
		make_dir("assembly_comparison");
		make_dir("assembly_comparison/saves");
		make_dir("assembly_comparison/initial_pics/");
		gp_t gp;
//		io::RCReaderWrapper<io::SingleRead> final_stream1(stream1_);
//		io::RCReaderWrapper<io::SingleRead> final_stream2(stream2_);
//		CompositeContigStream stream(final_stream1, final_stream2);
		CompositeContigStream stream(stream1_, stream2_);
		INFO("Constructing graph");
		ConstructGraph<gp_t::k_value, Graph>(gp.g, gp.index, stream);

		cout << "oppa" << endl;
		SimplifyGraph(gp.g);
		SplitGraph(gp);
		SaveGraph("assembly_comparison/saves/graph", gp);

		ColorHandler<Graph> coloring(gp.g);
		ColorGraph(gp, coloring);

		UntangledGraphConstructor<gp_t> ugp(gp, coloring, stream1_, stream2_);

//		INFO("Filling contig positions");
//		stream1_.reset();
//		FillPos(gp, stream1_, name1);
//		stream2_.reset();
//		FillPos(gp, stream2_, name2);

		SimplifyGraph(ugp.graph_);

		LengthIdGraphLabeler<Graph> labeler(ugp.graph_);

		ReliableSplitter<Graph> splitter(ugp.graph_, 30, 3000);
		ComponentSizeFilter<Graph> filter(ugp.graph_, 1000000000, 2);
		MapColorer<Graph, VertexId> vertex_colorer(ugp.untangled_coloring_.VertexColorMap());
		WriteComponents(ugp.graph_, splitter, filter, "breakpoint_graph",
				"assembly_comparison/initial_pics/breakpoint_graph.dot",
				ugp.untangled_coloring_.EdgeColorMap(), vertex_colorer, labeler);

		INFO("Removing unnecessary edges");
		DeleteVioletEdges(ugp.graph_, ugp.untangled_coloring_);

//		ReliableSplitter<Graph> splitter(gp.g, /*max_size*/100, /*edge_length_bound*/5000);
//		BreakPointsFilter<Graph> filter(gp.g, coloring, 3);
		INFO("Counting stats, outputting pictures");
		BPGraphStatCounter<Graph> counter(ugp.graph_, ugp.untangled_coloring_);
		counter.CountStats(labeler,
				MapColorer<Graph, VertexId>(ugp.untangled_coloring_.VertexColorMap()));

//		BPGraphStatCounter<Graph> counter(gp.g, coloring);
//		LengthIdGraphLabeler<Graph> labeler(gp.g);
//		counter.CountStats(labeler,
//				MapColorer<Graph, VertexId>(coloring.VertexColorMap()));
	}
private:
	DECL_LOGGER("AssemblyComparer")
	;
};

}
