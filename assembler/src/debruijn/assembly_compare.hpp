#pragma once

#include "standard.hpp"
#include "utils.hpp"
#include "omni_utils.hpp"
#include "graph_pack.hpp"
#include "graph_construction.hpp"

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
			vector<Range>& curr_ranges = crs[edge];
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
	typedef typename graph_pack::graph_t Graph;
	typedef typename Graph::EdgeId EdgeId;
	typedef map<EdgeId, vector<Range>> CoveredRanges;
	typedef map<EdgeId, vector<size_t>> BreakPoints;

	void FindBreakPoints(BreakPoints& bps, const CoveredRanges& crs1, const CoveredRanges& crs2) {
		map<EdgeId, set<size_t>> tmp;
		for (auto it = crs1.begin(); it != crs1.end(); ++it) {
			EdgeId edge = it->first;
			const vector<Range>& ranges = it->second;
			for (auto r_it = ranges.begin(); r_it != ranges.end(); ++r_it) {

			}
		}
	}

public:

	template<size_t k>
	void CompareAssemblies(ContigStream& ass1, ContigStream& ass2) {
		typedef graph_pack<ConjugateDeBruijnGraph, k> gp_t;
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
		FindBreakPoints(bps, crs1, crs2);
	}
};

}
