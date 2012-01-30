#pragma once

#include "standard.hpp"
#include "utils.hpp"
#include "omni_utils.hpp"
#include "graph_pack.hpp"
#include "graph_construction.hpp"
#include "graph_simplification.hpp"

namespace debruijn_graph {
typedef io::IReader<io::SingleRead> ContigStream;
typedef io::MultifileReader<io::SingleRead> CompositeContigStream;

template<class gp_t>
class CoveredRangesFinder {
	const gp_t& gp_;

	typedef typename graph_pack::graph_t Graph;
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
		}
		return answer;
	}

	void ProcessPath(const MappingPath<EdgeId>& path, CoveredRanges& crs) {
		for (size_t i = 0; i < path.size(); ++i) {
			auto mapping = path[i];
			EdgeId edge = mapping.first;
			const vector<Range>& curr_ranges = crs[edge];
			Range mapping_range = mapping.second.mapped_range;
			crs[edge] = ProcessRange(mapping_range, curr_ranges);
		}
	}

public:

	CoveredRangesFinder(const gp_t& gp) :
			gp_(gp) {

	}

	void FindCoveredRanges(CoveredRanges& crs, ContigStream& stream) {
		io::SingleRead read;
		NewExtendedSequenceMapper<gp_t::k_value, Graph> mapper(gp_.g, gp_.index,
				gp_.kmer_mapper);
		while (!stream.eof()) {
			stream >> read;
			ProcessPath(mapper.MapSequence(read.sequence()), crs);
		}
	}

};



template<class gp_t>
class AssemblyComparer {
public:
	enum edge_type {
		black = 0, red, blue, violet
	};
private:
	typedef typename graph_pack::graph_t Graph;
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
			while (breaks[i] != it->start_pos)
				++i;
			while (breaks[i] != it->end_pos) {
				coloring[i] += paint;
				++i;
			}
		}
	}

	pair<vector<size_t>, vector<edge_type>> CombineCoveredRanges(const vector<Range>& ranges1, const vector<Range>& ranges2) {
		set<size_t> tmp_breaks;
		AddBreaks(tmp_breaks, ranges1);
		AddBreaks(tmp_breaks, ranges2);
		vector<size_t> breaks(tmp_breaks.begin(), tmp_breaks.end());
		//breaks contain 0 and edge_length here!
		vector<edge_type> coloring;
		Paint(coloring, ranges1, breaks, edge_type::red);
		Paint(coloring, ranges1, breaks, edge_type::blue);
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
		}
	}

	void SplitEdge(const vector<size_t>& breaks, vector<edge_type>& colors, EdgeId e, Graph& g, map<EdgeId, string>& coloring) {
		VERIFY(breaks.size() + 1 == colors.size());
		vector<size_t> shifts;
		if (!breaks.empty) {
			shifts[0] = breaks[0];
			for (size_t i = 1; i < breaks.size(); ++i) {
				shifts[i] = breaks[i] - breaks[i-1];
			}
		}
		EdgeId curr_e = e;
		for (size_t i = 0; i < breaks.size(); ++i) {
			auto split_result = g.SplitEdge(curr_e, shifts[i]);
			coloring.insert(make_pair(split_result.first, colors[i]));
			curr_e = split_result.second;
		}
		coloring.insert(make_pair(curr_e, colors[breaks.size()]));
	}

	void SplitGraph(/*const */BreakPoints& bps, Graph& g, map<EdgeId, string>& coloring) {
		for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
			EdgeId e = *it;
			VERIFY(bps.find(e) != bps.end());
			SplitEdge(bps[e].first, bps[e].second, e, g, coloring);
		}
	}

public:

	template<size_t k>
	void CompareAssemblies(ContigStream& ass1, ContigStream& ass2) {
		gp_t gp;
		CompositeContigStream stream(ass1, ass2);
		ConstructGraph<k>(gp.g, gp.index, stream);
		CoveredRangesFinder<gp_t> crs_finder(gp);
		ass1.reset();
		ass2.reset();
		CoveredRanges crs1;
		crs_finder.FindCoveredRanges(crs1, ass1);
		CoveredRanges crs2;
		crs_finder.FindCoveredRanges(crs2, ass2);
		BreakPoints bps;
		FindBreakPoints(gp.g, bps, crs1, crs2);
		map<EdgeId, string> coloring;
		SplitGraph(bps, gp.g, coloring);
		RemoveBulges(gp.g, EmptyHandleF);

	}

};

}
