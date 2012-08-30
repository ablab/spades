#pragma once

#include "io/paired_read.hpp"
#include "seq_map.hpp"
#include "omni/omni_utils.hpp"
#include "omni/id_track_handler.hpp"
#include "logger/logger.hpp"
#include "omni/paired_info.hpp"
#include "xmath.h"
#include <boost/optional.hpp>
#include <iostream>
#include "sequence/sequence_tools.hpp"
#include "omni/splitters.hpp"

#include "runtime_k.hpp"

#include "kmer_map.hpp"
#include "new_debruijn.hpp"
//#include "common/io/paired_read.hpp"
namespace debruijn_graph {

template<class graph_pack, class read_type>
class MismatchShallNotPass {
private:
	typedef typename graph_pack::graph_t Graph;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef runtime_k::RtSeq Kmer;

	graph_pack &gp_;
	io::IReader<read_type>& stream_;
	double relative_threshold_;

	struct NuclCount {
		size_t counts_[4];
		NuclCount() {
      memset(counts_, 0, sizeof(counts_));
    }

		size_t &operator[](size_t nucl) {
			return counts_[nucl];
		}
	};

	typedef vector<NuclCount> MMInfo;

map<EdgeId, MMInfo> CountStatistics() {
	map<EdgeId, MMInfo> statistics;
	for(auto it = gp_.g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
		statistics[*it] = MMInfo(gp_.g.length(*it) + gp_.g.k());
	}
	stream_.reset();
	NewExtendedSequenceMapper<Graph> sm(gp_.g, gp_.index, gp_.kmer_mapper, gp_.g.k() + 1);
	while(!stream_.eof()) {
		read_type read;
		stream_ >> read;
		Sequence s_read = read.sequence();
		omnigraph::MappingPath<EdgeId> path = sm.MapSequence(s_read);
		if(path.size() == 1 && path[0].second.initial_range.size() == path[0].second.mapped_range.size()) {
			Range initial_range = path[0].second.initial_range;
			Range mapped_range = path[0].second.mapped_range;
			Sequence s_edge = gp_.g.EdgeNucls(path[0].first);
			size_t len = initial_range.size() + gp_.g.k();
			size_t cnt = 0;
			for(size_t i = 0; i < len; i++) {
				if(s_read[initial_range.start_pos + i] != s_edge[mapped_range.start_pos + i]) {
					cnt++;
				}
			}
			if(cnt <= 1) {
				MMInfo &info = statistics[path[0].first];
				for(size_t i = 0; i < len; i++) {
					size_t nucl_code = s_read[initial_range.start_pos + i];
					info[mapped_range.start_pos + i][nucl_code]++;
				}
			}
		}
	}
	return statistics;
}

bool CheckMapper(const Sequence &from, const Sequence &to) {
	return gp_.kmer_mapper.CheckCanRemap(from, to);
}

EdgeId CorrectNucl(EdgeId edge, size_t position, char nucl) {
	VERIFY(position >= gp_.g.k());
	if(position + 1 < gp_.g.length(edge)) {
		edge = gp_.g.SplitEdge(edge, position + 1).first;
	}
	EdgeId mismatch = edge;
	if(position > gp_.g.k()) {
		auto tmp = gp_.g.SplitEdge(edge, position - gp_.g.k());
		edge = tmp.first;
		mismatch = tmp.second;
	}
	Sequence s_mm = gp_.g.EdgeNucls(mismatch);
	VERIFY(nucl != s_mm[gp_.g.k()]);
	Sequence correct = s_mm.Subseq(0, gp_.g.k()) + Sequence(string(1, nucl)) + s_mm.Subseq(gp_.g.k() + 1, gp_.g.k() * 2 + 1);
	if(!CheckMapper(s_mm, correct)) {
		return edge;
	}
	EdgeId correct_edge = gp_.g.AddEdge(gp_.g.EdgeStart(mismatch), gp_.g.EdgeEnd(mismatch), correct);
	if(position > gp_.g.k()) {
		gp_.g.GlueEdges(mismatch, correct_edge);
		return edge;
	} else {
		return gp_.g.GlueEdges(mismatch, correct_edge);
	}
}

void CorrectNucls(EdgeId edge, vector<pair<size_t, char>> mismatches) {
	for(auto it = mismatches.rbegin(); it != mismatches.rend(); ++it) {
		edge = CorrectNucl(edge, it->first, it->second);
	}
	Compressor<Graph>(gp_.g).CompressVertex(gp_.g.EdgeEnd(edge));
}

vector<pair<size_t, char>> FindMismatches(EdgeId edge, MMInfo &statistics) {
	vector<pair<size_t, char>> to_correct;
	Sequence s_edge = gp_.g.EdgeNucls(edge);
	for(size_t i = gp_.g.k(); i < gp_.g.length(edge); i++) {
		size_t cur_best = 0;
		for(size_t j = 1; j < 4; j++) {
			if(statistics[i][j] > statistics[i][cur_best]) {
				cur_best = j;
			}
		}
		size_t nucl_code = s_edge[i];
		if(statistics[i][cur_best] > relative_threshold_ * statistics[i][nucl_code] + 1) {
			to_correct.push_back(make_pair(i, cur_best));
			i += gp_.g.k();
		}
	}
	return to_correct;
}

size_t CorrectEdge(EdgeId edge, MMInfo &statistics) {
	vector<pair<size_t, char>> to_correct = FindMismatches(edge, statistics);
	CorrectNucls(edge, to_correct);
	return to_correct.size();
}

public:
	MismatchShallNotPass(graph_pack &gp, io::IReader<read_type>& stream, double relative_threshold = 2) : gp_(gp), stream_(stream), relative_threshold_(relative_threshold) {
		VERIFY(relative_threshold >= 1);
	}

	size_t StopMismatch() {
		size_t res = 0;
		map<EdgeId, MMInfo> statistics = CountStatistics();
		for(auto it = gp_.g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
			if(!gp_.g.RelatedVertices(gp_.g.EdgeStart(*it), gp_.g.EdgeEnd(*it)) && statistics.find(*it) != statistics.end()) {
				res += CorrectEdge(*it, statistics[*it]);
			}
		}
		return res;
	}

	size_t StopAllMismatches(size_t max_iterations) {
		size_t res = 0;
		while(max_iterations > 0) {
			size_t last = StopMismatch();
			res += last;
			if(last == 0)
				break;
			max_iterations--;
		}
		return res;
	}
};

}
