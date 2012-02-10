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
		NewExtendedSequenceMapper<gp_t::k_value + 1, Graph> mapper(gp_.g, gp_.index,
				gp_.kmer_mapper);
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

	const map<EdgeId, string>& coloring_;
	size_t color_threshold_;
public:
	BreakPointsFilter(const Graph& graph, const map<EdgeId, string>& coloring, size_t color_threshold) :
			base(graph), coloring_(coloring), color_threshold_(color_threshold) {

	}

	bool MultiColored(const GraphComponent<Graph>& component) const {
		set<string> colors;
		for (auto it = component.e_begin(); it != component.e_end(); ++it) {
			auto color_it = coloring_.find(*it);
			VERIFY(color_it != coloring_.end());
			colors.insert(color_it->second);
		}
		return colors.size() >= color_threshold_;
	}

	/*virtual*/
	//todo change to set or GraphComponent and add useful protected methods
	bool Check(const vector<VertexId> &component_veritces) const {
		GraphComponent<Graph> component(this->graph(), component_veritces.begin(), component_veritces.end());
		return component.v_size() > 2 && MultiColored(component);
	}

};

template<class gp_t>
class AssemblyComparer {
public:
	enum edge_type {
		//don't change order!!!
		black = 0, red, blue, violet
	};

	static string color_str(edge_type color) {
		static string colors[] = {"black", "red", "blue", "purple"};
		return colors[(int)color];
	}

private:
	typedef typename gp_t::graph_t Graph;
	typedef typename Graph::EdgeId EdgeId;
	typedef map<EdgeId, vector<Range>> CoveredRanges;
	typedef map<EdgeId, pair<vector<size_t>, vector<edge_type>>> BreakPoints;

	void AddBreaks(set<size_t>& breaks, const vector<Range>& ranges) {
		for (auto it = ranges.begin(); it != ranges.end(); ++it) {
			breaks.insert(it->start_pos);
			breaks.insert(it->end_pos);
		}
	}

	void Paint(vector<edge_type>& coloring, const vector<Range>& ranges, const vector<size_t>& breaks, edge_type paint) {
		size_t i = 0;
		for (auto it = ranges.begin(); it != ranges.end(); ++it) {
			while (breaks[i] != it->start_pos) {
				++i;
				VERIFY(i < breaks.size());
			}
			while (breaks[i] != it->end_pos) {
				coloring[i] = (edge_type)((int)coloring[i] | (int) paint);
				++i;
				VERIFY(i < breaks.size());
			}
		}
	}

	pair<vector<size_t>, vector<edge_type>> CombineCoveredRanges(const vector<Range>& ranges1, const vector<Range>& ranges2) {
		set<size_t> tmp_breaks;
		AddBreaks(tmp_breaks, ranges1);
		AddBreaks(tmp_breaks, ranges2);
		vector<size_t> breaks(tmp_breaks.begin(), tmp_breaks.end());
		VERIFY(breaks.size() >= 2);
		//breaks contain 0 and edge_length here!
		vector<edge_type> coloring(breaks.size() - 1, edge_type::black);
		Paint(coloring, ranges1, breaks, edge_type::red);
		Paint(coloring, ranges2, breaks, edge_type::blue);
		//cleaning breaks from 0 and edge_length
		vector<size_t> final_breaks;
		for (size_t i = 1; i < breaks.size() - 1; ++i) {
			final_breaks.push_back(breaks[i]);
		}
		return make_pair(final_breaks, coloring);
	}

	void FindBreakPoints(const Graph& g, BreakPoints& bps, /*const */CoveredRanges& crs1, /*const */CoveredRanges& crs2) {
		for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
			EdgeId e = *it;
			bps[e] = CombineCoveredRanges(crs1[e], crs2[e]);
			VERIFY(bps[e].first.empty() || bps[e].first.back() < g.length(e));
		}
	}

	void SplitEdge(const vector<size_t>& breaks, const vector<edge_type>& colors, EdgeId e, Graph& g, map<EdgeId, string>& coloring) {
		VERIFY(breaks.size() + 1 == colors.size());
		vector<size_t> shifts(breaks.size());
		if (!breaks.empty()) {
			shifts[0] = breaks[0];
			for (size_t i = 1; i < breaks.size(); ++i) {
				shifts[i] = breaks[i] - breaks[i-1];
			}
		}
		EdgeId curr_e = e;
		for (size_t i = 0; i < breaks.size(); ++i) {
			auto split_result = g.SplitEdge(curr_e, shifts[i]);
			//todo doesn't work for paired graph!!!
			coloring.insert(make_pair(split_result.first, color_str(colors[i])));
			curr_e = split_result.second;
		}
		coloring.insert(make_pair(curr_e, color_str(colors[breaks.size()])));
	}

	void SplitGraph(/*const */BreakPoints& bps, Graph& g, map<EdgeId, string>& coloring) {
		set<EdgeId> initial_edges;
		for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
			initial_edges.insert(*it);
		}
		for (auto it = SmartSetIterator<Graph, EdgeId>(g, initial_edges.begin(), initial_edges.end()); !it.IsEnd(); ++it) {
			EdgeId e = *it;
			VERIFY(bps.find(e) != bps.end());
			VERIFY(bps[e].first.empty() || bps[e].first.back() < g.length(e));
			SplitEdge(bps[e].first, bps[e].second, e, g, coloring);
		}
	}

	void FillPos(gp_t& gp, ContigStream& stream, string stream_prefix) {
		stream.reset();
		io::SingleRead read;
		while (!stream.eof()) {
			stream >> read;
			FillEdgesPos(gp, read.sequence(), stream_prefix + read.name());
		}
	}

public:

	void CompareAssemblies(ContigStream& stream1, ContigStream& stream2
			, const string& name1, const string& name2) {
		gp_t gp;
		CompositeContigStream stream(stream1, stream2);
		ConstructGraph<gp_t::k_value, Graph>(gp.g, gp.index, stream);

		debruijn_config::simplification::bulge_remover br_config;
		br_config.max_bulge_length_coefficient = 4;
		br_config.max_coverage = 1000.;
		br_config.max_relative_coverage = 1.2;
		br_config.max_delta = 3;
		br_config.max_relative_delta = 0.1;
		RemoveBulges(gp.g, br_config);

		CoveredRangesFinder<gp_t> crs_finder(gp);
		CoveredRanges crs1;
		crs_finder.FindCoveredRanges(crs1, stream1);
		CoveredRanges crs2;
		crs_finder.FindCoveredRanges(crs2, stream2);
		BreakPoints bps;
		FindBreakPoints(gp.g, bps, crs1, crs2);
		map<EdgeId, string> coloring;

		//todo color after split!!! for dealing with conjugate graphs
		SplitGraph(bps, gp.g, coloring);

		FillPos(gp, stream1, name1);
		FillPos(gp, stream2, name2);

		EdgePosGraphLabeler<Graph> pos_labeler(gp.g, gp.edge_pos);
		StrGraphLabeler<Graph> str_labeler(gp.g);
		CompositeLabeler<Graph> labeler(pos_labeler, str_labeler);

		ReliableSplitter<Graph> splitter(gp.g, /*max_size*/100, /*edge_length_bound*/5000);
		BreakPointsFilter<Graph> filter(gp.g, coloring, 3);
		make_dir("assembly_comparison");
		WriteComponents(gp.g, splitter, filter,
				"breakpoint_graph", "assembly_comparison/breakpoint_graph.dot",
				coloring, labeler);
	}
};

template<class gp_t>
class ComponentClassificator {
private:
	typedef typename gp_t::graph_t Graph;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef typename AssemblyComparer<gp_t>::edge_type edge_type;

public:
	enum component_type {
		error = 0,
		single_red,
		single_blue,
		simple_bulge,
		simple_misassembly,
		complex_misassembly,
		size
	};

private:
	const gp_t &gp_;
	const size_t bulge_length_;

public:
	ComponentClassificator(const gp_t &gp, size_t bulge_length) :
			gp_(gp), bulge_length_(bulge_length) {
	}

	ComponentClassificator(const gp_t &gp) :
			gp_(gp), bulge_length_(gp.g.k() * 5) {
	}

	edge_type GetColour(EdgeId edge) {
		return edge_type::red;
	}

	bool CheckSimpleMisassembly(vector<size_t> sources, vector<size_t> sinks) {
		if(sources.size() != 2 || sinks.size() != 2)
			return false;
		for(size_t i = 0; i < sources.size(); i++)
			for(size_t j = 0; j < sinks.size(); j++) {
				if(gp_.g.GetEdgesBetween(sources[i], sinks[j]).size() != 1) {
					return false;
				}
			}
		for(size_t i = 0; i < 2; i++) {
			if(gp_.g.GetEdgesBetween(sources[i], sources[1 - i]).size() != 0)
				return false;
			if(gp_.g.GetEdgesBetween(sinks[i], sinks[1 - i]).size() != 0)
				return false;
		}
		for(size_t i = 0; i < 2; i++) {
			if(GetColour(gp_.g.GetEdgesBetween(sources[i], sinks[0])) == GetColour(gp_.g.GetEdgesBetween(sources[i], sinks[1])))
				return false;
			if(GetColour(gp_.g.GetEdgesBetween(sources[0], sinks[i])) == GetColour(gp_.g.GetEdgesBetween(sources[1], sinks[i])))
				return false;
		}
		return true;
	}

	component_type GetComponentType(const vector<VertexId> &component) {
		if(component.size() < 2)
			return component_type::error;
		if(component.size() == 2) {
			vector<EdgeId> edges01 = gp_.g.GetEdgesBetween(component[0], component[1]);
			vector<EdgeId> edges10 = gp_.g.GetEdgesBetween(component[0], component[1]);
			if(edges01.size() > 0 && edges10.size() > 0) {
				return component_type::complex_misassembly;
			}
			vector<EdgeId> edges;
			edges.insert(edges.end(), edges01.begin(), edges01.end());
			edges.insert(edges.end(), edges10.begin(), edges10.end());
			if(edges.size() == 1) {
				if(GetColour(edges[0]) == edge_type::red) {
					return component_type::single_red;
				} else {
					return component_type::single_blue;
				}
			}
			if(edges.size() == 2 && gp_.g.length(edges[0]) < bulge_length_ && gp_.g.length(edges[1]) < bulge_length_) {
				return component_type::simple_bulge;
			}
			return component_type::complex_misassembly;
		}
		if(component.size() == 4) {
			for(size_t i = 0; i < 4; i++)
				for(size_t j = i + 1; j < 4; j++) {
					vector<size_t> sources;
					sources.push_back(i);
					sources.push_back(j);
					vector<size_t> sinks;
					for(int k = 0; k < 4; k++) {
						if(k != i && k != j)
							sinks.push_back(k);
					}
					if(CheckSimpleMisassembly(sources, sinks)) {
						return component_type::simple_misassembly;
					}
				}
		}
		return component_type::complex_misassembly;
	}
};

template<class gp_t>
class BreakPointGraphStatistics : public GraphComponentFilter<Graph> {
private:
	typedef typename gp_t::graph_t Graph;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef typename AssemblyComparer<gp_t>::edge_type edge_type;
	typedef typename ComponentClassificator<gp_t>::component_type component_type;

	const gp_t &gp_;
	ComponentClassificator<gp_t> cc_;
	vector<size_t> ct_counter;
	bool ready_;

public:
	BreakPointGraphStatistics(const gp_t &gp) : gp_(gp), cc_(gp), ct_counter(component_type::size), ready_(false) {
	}

	bool Check(vector<VertexId> component) {
		component_type t = cc_.GetComponentType(component);
		ct_counter[t]++;
		return t == component_type::complex_misassembly;
	}

	void CountStats(const map<EdgeId, string> &coloring, const GraphLabeler<Graph>& labeler) {
		make_dir("assembly_comparison");
		WriteComponents(gp_.g, LongEdgesExclusiveSplitter<Graph>(gp_.g, 1000000000), *this,
			"breakpoint_graph", "assembly_comparison/breakpoint_graph.dot",
			coloring, labeler);
		ready_ = true;
	}

	size_t GetComponentNumber(component_type t) const {
		return ct_counter[t];
	}
};

}
