/*
 * genome_consistance_checker.hpp
 *
 *  Created on: Oct 24, 2013
 *      Author: anton
 */

#pragma once

#include "omni/edges_position_handler.hpp"
#include "omni/mapping_path.hpp"

namespace debruijn_graph {

template<class Graph>
class GenomeConsistenceChecker {
private:
	using omnigraph::MappingRange;
	const Graph &graph_;
	Sequence genome_;
	omnigraph::EdgesPositionHandler<Graph> genome_mapping_;
	size_t max_gap_;
	double relative_max_gap_;

	bool consequent(const MappingRange &mr1, const MappingRange &mr2, size_t gap) const {
		if (mr2.mapped_range.start_pos == 0 && mr1.mapped_range.end_pos == genome_.size())
			return true;
		if (mr1.initial_range.end_pos > mr2.initial_range.start_pos + gap)
			return false;
		if (mr1.initial_range.end_pos + gap > mr2.initial_range.start_pos)
			return false;
		if (mr1.mapped_range.end_pos > mr2.mapped_range.start_pos + gap)
			return false;
		if (mr1.mapped_range.end_pos + gap > mr2.mapped_range.start_pos)
			return false;
		return true;
	}

	set<MappingRange> FillPositionGaps(const set<MappingRange> &info, size_t gap) const {
		set<MappingRange> result;
		auto cur = info.begin();
		while(cur != info.end()) {
			MappingRange new_range = *cur;
			++cur;
			while(cur != info.end() && consequent(new_range, *cur, gap)) {
				new_range = new_range.Merge(*cur);
				++cur;
			}
			result.insert(new_range);
		}
		return result;
	}

	void Merge(set<MappingRange> &ranges, set<MappingRange> &to_merge, int shift) {
		for(set<MappingRange>::iterator it = to_merge.begin(); it != to_merge.end(); ++it) {
			ranges.insert(genome_mapping_.EraseAndExtract(ranges, it->Shift(shift)));
		}
	}

	bool IsConsistentWithGenomeStrand(const vector<EdgeId> &path, const string &strand) const {
		size_t len = graph_.length(path[0]);
		for (size_t i = 1; i < path.size(); i++) {
			Merge(res, genome_mapping_.GetEdgePositions(path[i], strand));
			len += graph_.length(path[i]);
		}
		FillPositionGaps(res, len);
		if (res.size() > 0) {
			for (size_t i = 0; i < res.size(); i++) {
				size_t m_len = res[i].initial_range.size();
				if (abs(int(res[i].initial_range.size()) - int(len)) < max(1.0 * max_gap_, 0.07 * len))
					return true;
			}
		}
		return false;
	}

public:
	template<class GraphPack>
	GenomeConsistenceChecker(const GraphPack gp, size_t max_gap, double relative_max_gap /*= 0.2*/) :
			graph_(gp.g), genome_(gp.genome), genome_mapping_(gp.g), max_gap_(max_gap), relative_max_gap_(relative_max_gap) {
        FillPos(gp, gp.genome, "0");
        FillPos(gp, !gp.genome, "1");
	}

	bool IsConsistentWithGenome(vector<EdgeId> path) const {
		if (path.size() == 0)
			return false;
		for (size_t i = 0; i + 1 < path.size(); i++) {
			if (graph_.EdgeStart(path[i + 1]) != graph_.EdgeEnd(path[i]))
				return false;
		}
		return IsConsistentWithGenomeStrand(path, "0") || IsConsistentWithGenomeStrand(path, "1");
	}
};
}
