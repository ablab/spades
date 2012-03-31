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

template<size_t k, class Graph>
void ConstructGraph(Graph& g, EdgeIndex<k + 1, Graph>& index,
		io::IReader<io::SingleRead>& stream1, io::IReader<io::SingleRead>& stream2) {
	io::MultifileReader<io::SingleRead> composite_reader(stream1, stream2);
	ConstructGraph<k, Graph>(g, index, composite_reader);
}

template<class Graph>
Sequence MergeSequences(const Graph& g,
		const vector<typename Graph::EdgeId>& continuous_path) {
	vector<Sequence> path_sequences;
	path_sequences.push_back(g.EdgeNucls(continuous_path[0]));
	for (size_t i = 1; i < continuous_path.size(); ++i) {
		VERIFY(g.EdgeEnd(continuous_path[i-1]) == g.EdgeStart(continuous_path[i]));
		path_sequences.push_back(g.EdgeNucls(continuous_path[i]));
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

	void AlignAndProject(const Sequence& tip_seq, const Sequence& alt_seq,
			bool outgoing_tip) {
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

		INFO("Remapping " << aligned_tip.size() << " kmers of aligned_tip to aligned_alt");
		gp_.kmer_mapper.RemapKmers(aligned_tip, aligned_alt);
	}

	void AlignAndProject(const AbstractConjugateGraph<typename Graph::DataMaster>& graph
			, const Sequence& tip_seq, const Sequence& alt_seq,
			bool outgoing_tip) {
		AlignAndProject(tip_seq, alt_seq, outgoing_tip);
		AlignAndProject(!tip_seq, !alt_seq, !outgoing_tip);
	}

	void AlignAndProject(const AbstractNonconjugateGraph<typename Graph::DataMaster>& graph
			, const Sequence& tip_seq, const Sequence& alt_seq,
			bool outgoing_tip) {
		AlignAndProject(tip_seq, alt_seq, outgoing_tip);
	}

public:
	TipsProjector(gpt& gp) :
			gp_(gp), unique_path_finder_(gp.g) {

	}

	void ProjectTip(EdgeId tip) {
		INFO("Trying to project tip " << gp_.g.str(tip));
		bool outgoing_tip = gp_.g.IsDeadEnd(gp_.g.EdgeEnd(tip));
		Sequence tip_seq = gp_.g.EdgeNucls(tip);
		vector<EdgeId> alt_path = UniqueAlternativePath(tip, outgoing_tip);
		if (alt_path.empty()) {
			WARN("Failed to find unique alt path for tip " << gp_.g.str(tip) <<". Wasn't projected!!!");
		} else {
			Sequence alt_seq = MergeSequences(gp_.g, alt_path);
			if (tip_seq.size() > alt_seq.size()) {
				WARN("Can't fully project tip " << gp_.g.str(tip) << " with seq length " << tip_seq.size()
						<< " because alt path length is " << alt_seq.size() << ". Trying to project partially");
			}
			AlignAndProject(gp_.g, tip_seq, alt_seq, outgoing_tip);
			INFO("Tip projected");
		}
	}
private:
	DECL_LOGGER("TipsProjector")
	;
};

//beware!!! this class clears quality!!!
class ModifyingWrapper: public io::DelegatingReaderWrapper<io::SingleRead> {
protected:

	ModifyingWrapper(io::IReader<io::SingleRead>& reader) :
			io::DelegatingReaderWrapper<io::SingleRead>(reader) {
	}

	virtual Sequence Modify(const Sequence& read) const = 0;

public:
	ModifyingWrapper& operator>>(io::SingleRead& read) {
		this->reader() >> read;
		read = io::SingleRead(read.name(), Modify(read.sequence()).str());
		return *this;
	}
};

template<class gp_t>
class ContigRefiner: public ModifyingWrapper {
	typedef ModifyingWrapper base;
	typedef typename gp_t::graph_t Graph;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	static const size_t k = gp_t::k_value;
	const Graph& graph_;
	NewExtendedSequenceMapper<k + 1, Graph> mapper_;

	//todo seems that we don't need optional here any more
	boost::optional<vector<EdgeId>> TryCloseGap(VertexId v1, VertexId v2) const {
		PathStorageCallback<Graph> path_store(graph_);
		//todo reduce value after investigation
		PathProcessor<Graph> path_processor(graph_, 0, 50, v1, v2, path_store);
		path_processor.Process();

		if (path_store.count() == 0) {
			VERIFY_MSG(false, "Failed to find closing path");
//			return boost::none;
		} else if (path_store.count() == 1) {
			return boost::optional<vector<EdgeId>>(*path_store.paths().begin());
		} else {
//			VERIFY_MSG(false, "Several closing paths found");
			WARN("Several closing paths found");
			return boost::optional<vector<EdgeId>>(*path_store.paths().begin());
		}
	}

	vector<EdgeId> FixPath(const vector<EdgeId>& edges) const {
		vector<EdgeId> answer;
		if (edges.empty()) {
			WARN("Mapping path was empty");
			return vector<EdgeId>();
		}
//		VERIFY(edges.size() > 0);
		answer.push_back(edges[0]);
		for (size_t i = 1; i < edges.size(); ++i) {
			VertexId v1 = graph_.EdgeEnd(edges[i - 1]);
			VertexId v2 = graph_.EdgeStart(edges[i]);
			if (v1 == v2) {
				answer.push_back(edges[i]);
			} else {
				INFO(
						"Trying to close gap between v1=" << graph_.int_id(v1)
								<< " and v2=" << graph_.int_id(v2));
				boost::optional<vector<EdgeId>> closure = TryCloseGap(v1, v2);
				if (closure) {
					INFO("Gap closed");
					INFO("Cumulative closure length is " << CummulativeLength(graph_, *closure));
					answer.insert(answer.end(), closure->begin(),
							closure->end());
					answer.push_back(edges[i]);
				} else {
					WARN(
							"Failed to close gap between v1="
									<< graph_.int_id(v1) << " (conjugate "
									<< graph_.int_id(graph_.conjugate(v1))
									<< ") and v2="
									<< graph_.int_id(v2) << " (conjugate "
									<< graph_.int_id(graph_.conjugate(v2))
									<< ")");
					make_dir("assembly_compare/tmp");
					WriteComponentsAroundEdge(graph_, graph_.IncomingEdges(v1).front()
							, "assembly_compare/tmp/failed_close_gap_from.dot"
							, *DefaultColorer(graph_), LengthIdGraphLabeler<Graph>(graph_));
					WriteComponentsAroundEdge(graph_, graph_.OutgoingEdges(v2).front()
							, "assembly_compare/tmp/failed_close_gap_to.dot"
							, *DefaultColorer(graph_), LengthIdGraphLabeler<Graph>(graph_));

					VERIFY(false);
				}
			}
//			VERIFY(graph_.EdgeStart(path[i]) == graph_.EdgeEnd(path[i - 1]));
		}
		return answer;
	}

	Path<EdgeId> FixPath(const Path<EdgeId>& path) const {
		return Path<EdgeId>(FixPath(path.sequence()), path.start_pos(),
				path.end_pos());
	}

protected:

	Sequence Modify(const Sequence& s) const {
//		if(s < !s)
//			return !Refine(!s);
		Path<EdgeId> path = FixPath(mapper_.MapSequence(s).simple_path());

		if (path.size() == 0) {
			WARN("For sequence of length " << s.size() << " returning empty sequence");
			return Sequence();
		}
//		DEBUG("Mapped sequence to path " << graph_.str(path.sequence()));

		Sequence path_sequence = MergeSequences(graph_, path.sequence());
		size_t start = path.start_pos();
		size_t end = path_sequence.size()
				- graph_.length(path[path.size() - 1]) + path.end_pos();
		//todo output levenshtein somewhere there!!
		Sequence answer = path_sequence.Subseq(start, end);
		if (answer != s) {
			if (answer.size() < 1000) {
				TRACE("Initial sequence modified, edit distance= " << EditDistance(answer, s));
			} else {
				TRACE("Sequence too large, won't count edit distance");
			}
		}
//		else {
//			TRACE("Initial sequence unmodified!");
//		}
		return answer;
	}

public:
	ContigRefiner(io::IReader<io::SingleRead>& reader, const gp_t& gp) :
			base(reader), graph_(gp.g), mapper_(gp.g, gp.index, gp.kmer_mapper) {
	}

//	/* virtual */
//	ContigRefiner& operator>>(io::SingleRead& read) {
//		this->reader() >> read;
//		Sequence s = Refine(read.sequence());
//		read.SetSequence(s.str().c_str());
//		string quality(s.size(), (char)33);
//		read.SetQuality(quality.c_str());
//		return *this;
//	}
private:
	DECL_LOGGER("ContigRefiner");
};

template<class gp_t>
void ConstructGPForRefinement(gp_t& gp,
		io::IReader<io::SingleRead>& raw_stream_1,
		io::IReader<io::SingleRead>& raw_stream_2, size_t delta) {
	typedef typename gp_t::graph_t Graph;
	INFO("Constructing graph pack for refinement");

	io::RCReaderWrapper < io::SingleRead > stream_1(raw_stream_1);
	io::RCReaderWrapper < io::SingleRead > stream_2(raw_stream_2);
	stream_1.reset();
	stream_2.reset();

	ConstructGraph<gp_t::k_value>(gp.g, gp.index, stream_1, stream_2);

//	make_dir("bp_graph_test/tmp/");
//	LengthIdGraphLabeler<Graph> labeler(gp.g);
//	WriteToDotFile(gp.g, labeler, "bp_graph_test/tmp/before_refine.dot");

	//todo configure!!!
	debruijn_config::simplification::bulge_remover br_config;
	br_config.max_bulge_length_coefficient = 3;
	br_config.max_coverage = 1000.;
	br_config.max_relative_coverage = 1.2;
	br_config.max_delta = delta;
	br_config.max_relative_delta = 0.1;

	INFO("Removing bulges");
	RemoveBulges(gp.g, br_config);

	INFO("Remapped " << gp.kmer_mapper.size() << " k-mers");

	TipsProjector<gp_t> tip_projector(gp);
	boost::function<void(EdgeId)> projecting_callback = boost::bind(
			&TipsProjector<gp_t>::ProjectTip, &tip_projector, _1);
	debruijn_config::simplification::tip_clipper tc_config;
	tc_config.max_coverage = 1000.;
	tc_config.max_relative_coverage = 1.1;
	tc_config.max_tip_length_coefficient = 2.;

	INFO("Clipping tips with projection");
	ClipTips(gp.g, tc_config, /*read_length*/100, projecting_callback);

	INFO("Remapped " << gp.kmer_mapper.size() << " k-mers");

//	WriteToDotFile(gp.g, labeler, "bp_graph_test/tmp/after_refine.dot");
}

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

	/*virtual*/ void HandleMerge(const vector<EdgeId>& old_edges, EdgeId new_edge) {
		VERIFY(old_edges.size() > 0);
		auto color = Color(old_edges.front());
		for (auto it = old_edges.begin(); it != old_edges.end(); ++it) {
			VERIFY(color == Color(*it));
//			PaintEdge(new_edge, Color(*it));
		}
		Paint(new_edge, color);
	}
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
		GraphComponent < Graph
				> component(this->graph(), component_veritces.begin(),
						component_veritces.end());
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

	size_t to_draw_;
	ComponentClassifier<Graph> cc_;

public:
	ComponentTypeFilter(const Graph &g, size_t to_draw,
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
auto_ptr<GraphColorer<Graph>> ConstructColorer(
		const ColorHandler<Graph>& coloring) {
	return auto_ptr < GraphColorer
			< Graph
					>> (new CompositeGraphColorer<Graph>(
							new MapColorer<Graph, typename Graph::VertexId>(
									coloring.VertexColorMap())
							, new MapColorer<Graph, typename Graph::EdgeId>(
									coloring.EdgeColorMap())));
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

	void PrintComponents(size_t c_type, const GraphLabeler<Graph>& labeler,
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

class RCSplittingStream: public io::DelegatingReaderWrapper<io::SingleRead> {
private:
	io::SplittingWrapper filtered_reader_;
	io::RCReaderWrapper<io::SingleRead> stream_;
public:
	RCSplittingStream(ContigStream &base_stream) :
			filtered_reader_(base_stream), stream_(filtered_reader_) {
		Init(stream_);
	}
};

template<class Graph>
struct bp_graph_pack {
	typedef Graph graph_t;
	typedef string contig_id_t;
	typedef typename Graph::EdgeId EdgeId;
	Graph g;
	IdTrackHandler<Graph> int_ids;
	ColorHandler<Graph> coloring;
	map<contig_id_t, vector<EdgeId>> red_paths;
	map<contig_id_t, vector<EdgeId>> blue_paths;
	EdgesPositionHandler<Graph> edge_pos;

	bp_graph_pack(size_t k) :
			g(k), int_ids(g), coloring(g), edge_pos(g) {

	}
};

//todo finish later
//template<class gp_t1, class gp_t2>
//void ConvertToBPGraphPack(const gp_t1& gp
//		, const ColorHandler<typename gp_t1::graph_t>& coloring
//		, gp_t2& bp_gp) {
//	string tmp_dir = "/home/snurk/tmp/";
//	string filename = tmp_dir + "tmp";
//	make_dir(tmp_dir);
//	PrintGraphPack(filename, gp);
//	typename ScannerTraits<typename gp_t2::graph_t>::Scanner scanner(bp_gp.g,
//				bp_gp.int_ids);
//	ScanBasicGraph(filename, scanner);
//	scanner.loadPositions(filename, bp_gp.edge_pos);
//	//
//}

template<class gp_t>
class UntangledGraphContigMapper {
	typedef typename gp_t::graph_t Graph;
	typedef typename Graph::EdgeId EdgeId;
	const bp_graph_pack<Graph>& bp_gp_;

	MappingRange TrivialRange(EdgeId e, size_t& offset) const {
		size_t l = bp_gp_.g.length(e);
		offset += l;
		return MappingRange(Range(offset - l, offset), Range(0, 1));
	}

	MappingPath<EdgeId> TrivialMappingPath(const vector<EdgeId>& edges) const {
		vector<MappingRange> ranges;
		size_t offset = 0;
		for (auto it = edges.begin(); it != edges.end(); ++it) {
			ranges.push_back(TrivialRange(*it, offset));
		}
		return MappingPath<EdgeId>(edges, ranges);
	}

public:
	UntangledGraphContigMapper(const bp_graph_pack<Graph>& bp_gp) :
			bp_gp_(bp_gp) {

	}

	MappingPath<EdgeId> MapRead(const io::SingleRead &read) const {
		auto it = bp_gp_.red_paths.find(read.name());
		if (it != bp_gp_.red_paths.end()) {
			return TrivialMappingPath(it->second);
		}
		it = bp_gp_.blue_paths.find(read.name());
		if (it != bp_gp_.blue_paths.end()) {
			return TrivialMappingPath(it->second);
		}
		VERIFY(false);
		return MappingPath<EdgeId>();
	}

};

template<class gp_t>
class UntangledGraphConstructor {
private:
	typedef typename gp_t::graph_t Graph;
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	static const size_t k = gp_t::k_value;

	const gp_t& old_gp_;
	const ColorHandler<Graph>& old_coloring_;
	bp_graph_pack<Graph>& new_gp_;
	restricted::map<EdgeId, EdgeId> purple_edge_mapping_;
	restricted::map<VertexId, VertexId> vertex_mapping_;
	//todo draw in different color!
	restricted::set<VertexId> artificial_vertices_;
	set<string> processed_contigs_;

	//todo test that!!!
	string ConjugateContigId(const string& contig_id) {
		string answer;
		if (contig_id.substr(contig_id.size() - 3, 3) == "_RC")
			answer = contig_id.substr(0, contig_id.size() - 3);
		else
			answer = contig_id + "_RC";
		DEBUG("Conjugate to " << contig_id << " is " << answer);
		return answer;
	}

	void AddToProcessed(const string& contig_id) {
		processed_contigs_.insert(contig_id);
		processed_contigs_.insert(ConjugateContigId(contig_id));
	}

	VertexId GetStartVertex(const Path<EdgeId> &path, size_t i) {
		if (i != 0 || path.start_pos() == 0)
			return vertex_mapping_[old_gp_.g.EdgeStart(path[i])];
		else {
			//todo discuss with Anton!!!
			VertexId art_v = new_gp_.g.AddVertex();
			WARN("Art vertex added")
//			VERIFY(false);
			artificial_vertices_.insert(art_v);
			return art_v;
		}
	}

	VertexId GetEndVertex(const Path<EdgeId> &path, size_t i) {
		if (i != path.size() - 1 || path.end_pos() == old_gp_.g.length(path[i]))
			return vertex_mapping_[old_gp_.g.EdgeEnd(path[i])];
		else {
			//todo discuss with Anton!!!
			VertexId art_v = new_gp_.g.AddVertex();
			WARN("Art vertex added")
//			VERIFY(false);
			artificial_vertices_.insert(art_v);
			return art_v;
		}
	}

	void Untangle(ContigStream& stream, edge_type color) {
		io::SingleRead read;
		stream.reset();
		set<string> processed;
		while (!stream.eof()) {
			stream >> read;
			//todo can look at new_gp_.*_paths keys
			if (processed.count(read.name()) > 0)
				continue;
			processed.insert(read.name());
			processed.insert(ConjugateContigId(read.name()));

			Untangle(read.sequence(), read.name(), color);
		}
	}

	void Untangle(const Sequence& contig, const string& name, edge_type color) {
		VERIFY(color == edge_type::red || color == edge_type::blue);
		DEBUG("Untangling contig " << name);
		NewExtendedSequenceMapper<k + 1, Graph> mapper(old_gp_.g, old_gp_.index,
				old_gp_.kmer_mapper);
		Path<EdgeId> path = mapper.MapSequence(contig).simple_path();
		vector<EdgeId> new_path;
		DEBUG("Mapped contig" << name);
		for (size_t i = 0; i < path.size(); i++) {
			EdgeId next;
			if (old_coloring_.Color(path[i]) != edge_type::violet) {
				DEBUG("Next edge is not purple");
				size_t j = i;
				vector<EdgeId> to_glue;
				while (j < path.size()
						&& old_coloring_.Color(path[j]) != edge_type::violet) {
					to_glue.push_back(path[j]);
					j++;
				}
				Sequence new_edge_sequence = MergeSequences(old_gp_.g, to_glue);
				next = new_gp_.g.AddEdge(GetStartVertex(path, i),
						GetEndVertex(path, j - 1), new_edge_sequence);
				DEBUG("Added shortcut edge " << new_gp_.g.int_id(next) << " for path " << old_gp_.g.str(to_glue));
				i = j - 1;
			} else {
				DEBUG("Next edge is purple");
				next = purple_edge_mapping_[path[i]];
			}
			new_path.push_back(next);
			DEBUG("Coloring new edge and complement");
			PaintEdgeWithVertices(next, color);
		}
		if (color == edge_type::red) {
			VERIFY(new_gp_.red_paths.find(name) == new_gp_.red_paths.end());
			new_gp_.red_paths[name] = new_path;
			new_gp_.red_paths[ConjugateContigId(name)] = ConjugatePath(new_gp_.g, new_path);
		} else {
			VERIFY(new_gp_.blue_paths.find(name) == new_gp_.blue_paths.end());
			new_gp_.blue_paths[name] = new_path;
			new_gp_.blue_paths[name] = ConjugatePath(new_gp_.g, new_path);
		}
	}

	vector<EdgeId> ConjugatePath(const Graph& g, const vector<EdgeId> path) {
		vector<EdgeId> answer;
		for (int i = path.size() - 1; i >= 0; i--) {
			answer.push_back(g.conjugate(path[i]));
		}
		return answer;
	}

	template <class T>
	void ColorWithConjugate(T t, edge_type color) {
		new_gp_.coloring.Paint(t, color);
		new_gp_.coloring.Paint(new_gp_.g.conjugate(t), color);
	}

	void PaintEdgeWithVertices(EdgeId e, edge_type color) {
		DEBUG("Coloring edges " << new_gp_.g.int_id(e) << " and " << new_gp_.g.int_id(new_gp_.g.conjugate(e)));
		ColorWithConjugate(e, color);
		ColorWithConjugate(new_gp_.g.EdgeStart(e), color);
		ColorWithConjugate(new_gp_.g.EdgeEnd(e), color);
	}

public:
	UntangledGraphConstructor(const gp_t &old_gp,
			const ColorHandler<Graph> &old_coloring,
			bp_graph_pack<Graph>& new_gp, io::IReader<io::SingleRead> &stream1,
			io::IReader<io::SingleRead> &stream2) :
			old_gp_(old_gp), old_coloring_(old_coloring), new_gp_(new_gp) {
		const Graph& old_graph = old_gp.g;
		//adding vertices
		set<VertexId> processed_purple_v;
		for (auto it = old_graph.begin(); it != old_graph.end(); ++it) {
			if (processed_purple_v.count(*it) > 0)
				continue;
			processed_purple_v.insert(*it);
			processed_purple_v.insert(old_graph.conjugate(*it));
			vertex_mapping_[*it] = new_gp_.g.AddVertex();
			vertex_mapping_[old_graph.conjugate(*it)] = new_gp_.g.conjugate(vertex_mapping_[*it]);
			DEBUG("Adding purple vertex " << new_gp_.g.int_id(vertex_mapping_[*it])
					<< " corresponding to " << old_graph.int_id(*it) << " and conjugates")
		}


		set<EdgeId> processed_purple;
		//propagating purple color to new graph
		for (auto it = old_graph.SmartEdgeBegin(); !it.IsEnd(); ++it) {
			if (processed_purple.count(*it) > 0)
				continue;
			processed_purple.insert(*it);
			processed_purple.insert(old_graph.conjugate(*it));

			if (old_coloring.Color(*it) == edge_type::violet) {
				EdgeId new_edge = new_gp_.g.AddEdge(
						vertex_mapping_[old_graph.EdgeStart(*it)],
						vertex_mapping_[old_graph.EdgeEnd(*it)],
						old_graph.EdgeNucls(*it));
				DEBUG("Adding purple edge " << new_gp_.g.int_id(new_edge) << " corresponding to " << old_graph.int_id(*it) << " and conjugate")
				purple_edge_mapping_[*it] = new_edge;
				purple_edge_mapping_[old_graph.conjugate(*it)] = new_gp_.g.conjugate(new_edge);
				PaintEdgeWithVertices(new_edge, edge_type::violet);
			}
		}

		VERIFY(new_gp_.red_paths.empty());
		VERIFY(new_gp_.blue_paths.empty());

		Untangle(stream1, edge_type::red);
		Untangle(stream2, edge_type::blue);

		UntangledGraphContigMapper<bp_graph_pack<Graph>> contig_mapper(new_gp_);
		FillPos(new_gp_.g, contig_mapper, new_gp_.edge_pos, stream1);
		FillPos(new_gp_.g, contig_mapper, new_gp_.edge_pos, stream2);
	}
private:
	DECL_LOGGER("UntangledGraphConstructor");
};

template<class gp_t>
class AssemblyComparer {
private:
	typedef typename gp_t::graph_t Graph;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef map<EdgeId, vector<Range>> CoveredRanges;
	typedef map<EdgeId, vector<size_t>> BreakPoints;

//	io::IReader<io::SingleRead> &stream1_;
//	io::IReader<io::SingleRead> &stream2_;
	gp_t gp_;
	io::RCReaderWrapper<io::SingleRead> rc_stream1_;
	io::RCReaderWrapper<io::SingleRead> rc_stream2_;
	PrefixAddingReaderWrapper stream1_;
	PrefixAddingReaderWrapper stream2_;
	bool untangle_;

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

	void FindBreakPoints(BreakPoints& bps, /*const */
	CoveredRanges& crs1, /*const */CoveredRanges& crs2) {
		for (auto it = gp_.g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
			EdgeId e = *it;
			bps[e] = CombineCoveredRanges(crs1[e], crs2[e]);
			VERIFY(bps[e].empty() || bps[e].back() < gp_.g.length(e));
		}
	}

	void SplitEdge(const vector<size_t>& breaks, EdgeId e) {
		vector<size_t> shifts(breaks.size());
		if (!breaks.empty()) {
			shifts[0] = breaks[0];
			for (size_t i = 1; i < breaks.size(); ++i) {
				shifts[i] = breaks[i] - breaks[i - 1];
			}
		}
		EdgeId curr_e = e;
		for (size_t i = 0; i < breaks.size(); ++i) {
			auto split_result = gp_.g.SplitEdge(curr_e, shifts[i]);
			curr_e = split_result.second;
		}
	}

	void PrintCRS(const map<EdgeId, vector<Range>>& crs) {
		for (auto it = crs.begin(); it != crs.end(); ++it) {
			DEBUG(
					"For edge " << gp_.g.str(it->first) << " ranges "
							<< it->second);
		}
	}

	void SplitGraph() {
		INFO("Determining covered ranges");
		CoveredRangesFinder<gp_t> crs_finder(gp_);
		CoveredRanges crs1;
		crs_finder.FindCoveredRanges(crs1, stream1_);
		DEBUG("Printing covered ranges for stream 1");
		PrintCRS(crs1);
		CoveredRanges crs2;
		crs_finder.FindCoveredRanges(crs2, stream2_);
		DEBUG("Printing covered ranges for stream 2");
		PrintCRS(crs2);
		BreakPoints bps;
		INFO("Determining breakpoints");
		FindBreakPoints(bps, crs1, crs2);

		INFO("Splitting graph");
		SplitGraph(bps);
	}

	void SplitGraph(/*const */BreakPoints& bps) {
		set<EdgeId> initial_edges;
		for (auto it = gp_.g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
			initial_edges.insert(*it);
		}
		for (auto it = SmartSetIterator<Graph, EdgeId, typename Graph::Comparator>(gp_.g,
				initial_edges.begin(), initial_edges.end(), gp_.g.ReliableComparatorInstance()); !it.IsEnd();
				++it) {
			EdgeId e = *it;
			VERIFY(bps.find(e) != bps.end());
			VERIFY(bps[e].empty() || bps[e].back() < gp_.g.length(e));
			//todo temporary fix!!!
			if (e == gp_.g.conjugate(e))
				continue;
			SplitEdge(bps[e], e);
		}
	}

	void CompressGraph(Graph& g, ColorHandler<Graph>& coloring) {
		for (auto it = g.SmartVertexBegin(); !it.IsEnd(); ++it) {
			VertexId v = *it;
			if (g.CanCompressVertex(v)
					&& coloring.Color(g.GetUniqueOutgoingEdge(v))
							== coloring.Color(g.GetUniqueIncomingEdge(v))) {
				g.CompressVertex(v);
			}
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

	void ColorPath(const Path<EdgeId>& path, ColorHandler<Graph>& coloring,
			edge_type color) {
		for (size_t i = 0; i < path.size(); ++i) {
			coloring.Paint(path[i], color);
			coloring.Paint(gp_.g.EdgeStart(path[i]), color);
			coloring.Paint(gp_.g.EdgeEnd(path[i]), color);
		}
	}

	void ColorGraph(ColorHandler<Graph>& coloring) {
		INFO("Coloring graph");
		ColorGraph(coloring, stream1_, edge_type::red);
		ColorGraph(coloring, stream2_, edge_type::blue);
	}

	void ColorGraph(ColorHandler<Graph>& coloring, ContigStream& stream,
			edge_type color) {
		io::SingleRead read;
		stream.reset();
		NewExtendedSequenceMapper<gp_t::k_value + 1, Graph> mapper(gp_.g,
				gp_.index, gp_.kmer_mapper);
		while (!stream.eof()) {
			stream >> read;
			ColorPath(mapper.MapSequence(read.sequence()).simple_path(),
					coloring, color);
		}
	}

	void SimplifyGraph(Graph& g) {
		debruijn_config::simplification::bulge_remover br_config;
		br_config.max_bulge_length_coefficient = 50;
		br_config.max_coverage = 1000.;
		br_config.max_relative_coverage = 1.2;
		br_config.max_delta = 20;
		br_config.max_relative_delta = 0.1;
		INFO("Removing bulges");
		RemoveBulges(g, br_config);

//		debruijn_config::simplification::tip_clipper tc;
//		tc.max_coverage = 1000;
//		tc.max_relative_coverage = 1000;
//		tc.max_tip_length_coefficient = 6;
//		ClipTips(gp.g, tc, 10 * gp.g.k());
	}

	void SaveOldGraph(const string& path) {
		INFO("Saving graph to " << path);
		PrintGraphPack(path, gp_);
//		LengthIdGraphLabeler<Graph> labeler(gp_.g);
//		WriteToDotFile(gp_.g, labeler, path + ".dot");
	}

	template <class gp_t2>
	void UniversalSaveGP(gp_t2& gp, const string& filename) {
		typename PrinterTraits<Graph>::Printer printer(gp.g,
				gp.int_ids);
		INFO("Saving graph to " << filename);
		printer.saveGraph(filename);
		printer.saveEdgeSequences(filename);
		printer.savePositions(filename, gp.edge_pos);

		LengthIdGraphLabeler<Graph> labeler(gp.g);
		WriteToDotFile(gp.g, labeler, filename + ".dot");
	}

	void PrintColoredGraph(const Graph& g, const ColorHandler<Graph>& coloring,
			const EdgesPositionHandler<Graph>& pos,
			const string& output_filename) {
		ReliableSplitter<Graph> splitter(g, 30, 1000000);
		LengthIdGraphLabeler<Graph> basic_labeler(g);
		EdgePosGraphLabeler<Graph> pos_labeler(g, pos);

		CompositeLabeler<Graph> labeler(basic_labeler, pos_labeler);
		WriteComponents(g, splitter, output_filename,
				*ConstructColorer(coloring), labeler);
	}

	template <class gp_t2>
	void ProduceResults(gp_t2& gp, const ColorHandler<Graph>& coloring,
			const string& output_folder, bool detailed_output) {
		INFO("Removing unnecessary edges");
		DeleteVioletEdges(gp.g, coloring);

		if (detailed_output) {
			PrintColoredGraph(gp.g, coloring, gp.edge_pos,
					output_folder + "initial_pics/purple_removed.dot");
			UniversalSaveGP(gp, output_folder + "saves/purple_removed");
		}

//		ReliableSplitter<Graph> splitter(gp.g, /*max_size*/100, /*edge_length_bound*/5000);
//		BreakPointsFilter<Graph> filter(gp.g, coloring, 3);
		INFO("Counting stats, outputting pictures");
		BPGraphStatCounter<Graph> counter(gp.g, coloring, output_folder);
		LengthIdGraphLabeler<Graph> labeler(gp.g);
		counter.CountStats(labeler, detailed_output);
	}

public:

	AssemblyComparer(io::IReader<io::SingleRead> &stream1,
			io::IReader<io::SingleRead> &stream2, const string& name1,
			const string& name2, bool untangle = true) :
			rc_stream1_(stream1), rc_stream2_(stream2), stream1_(rc_stream1_,
					name1), stream2_(rc_stream2_, name2), untangle_(untangle) {
	}

	void CompareAssemblies(const string& output_folder, bool detailed_output =
			true) {
		//todo ???
		stream1_.reset();
		stream2_.reset();

		make_dir(output_folder);
		if (detailed_output) {
			make_dir(output_folder + "initial_pics/");
			make_dir(output_folder + "saves/");
			make_dir(output_folder + "purple_edges_pics/");
		}

		INFO("Constructing graph");
		INFO("K = " << gp_t::k_value);
		ConstructGraph<gp_t::k_value, Graph>(gp_.g, gp_.index, stream1_,
				stream2_);

		//TODO do we still need it?
//		SimplifyGraph(gp.g);

		if (detailed_output) {
			//saving for debug and investigation
			SaveOldGraph(output_folder + "saves/init_graph");
		}

		SplitGraph();

		if (detailed_output) {
			//saving for debug and investigation
			SaveOldGraph(output_folder + "saves/split_graph");
		}

		ColorHandler<Graph> coloring(gp_.g);
		ColorGraph(coloring);

		//situation in example 6 =)
		CompressGraph(gp_.g, coloring);

		INFO("Filling contig positions");
		stream1_.reset();
		FillPos<gp_t>(gp_, stream1_);
		stream2_.reset();
		FillPos<gp_t>(gp_, stream2_);

		if (detailed_output) {
			PrintColoredGraph(gp_.g, coloring, gp_.edge_pos,
					output_folder + "initial_pics/colored_split_graph.dot");
			SaveOldGraph(output_folder + "saves/tangled_graph");
		}

		if (untangle_) {
			INFO("Untangling graph");
			bp_graph_pack<typename gp_t::graph_t> untangled_gp(gp_t::k_value);
			UntangledGraphConstructor<gp_t> untangler(gp_, coloring, untangled_gp,
					stream1_, stream2_);
			//todo ???
			//		SimplifyGraph(untangled_gp.g);

			if (detailed_output) {
				PrintColoredGraph(untangled_gp.g, untangled_gp.coloring, untangled_gp.edge_pos,
						output_folder + "initial_pics/untangled_graph.dot");
				UniversalSaveGP(untangled_gp, output_folder + "saves/untangled_graph");
			}

			ProduceResults(untangled_gp, untangled_gp.coloring, output_folder,
					detailed_output);
		} else {
			ProduceResults(gp_, coloring, output_folder, detailed_output);
		}
	}

private:
	DECL_LOGGER("AssemblyComparer")
	;
};

}
