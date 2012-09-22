#pragma once

#include "io/paired_read.hpp"
#include "omni/omni_utils.hpp"
#include "omni/id_track_handler.hpp"
#include "logger/logger.hpp"
#include "omni/paired_info.hpp"
#include "xmath.h"
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
	double relative_threshold_;

	struct NuclCount {
		size_t counts_[4];
		NuclCount() {
      memset(counts_, 0, sizeof(counts_));
    }

		size_t &operator[](size_t nucl) {
			return counts_[nucl];
		}

		NuclCount &operator+=(const NuclCount &other) {
			counts_[0] += other.counts_[0];
			counts_[1] += other.counts_[1];
			counts_[2] += other.counts_[2];
			counts_[3] += other.counts_[3];
			return *this;
		}
	};

	typedef vector<NuclCount> MMInfo;

void CountStatistics(io::IReader<read_type>& stream, map<EdgeId, MMInfo> &statistics) {
	for(auto it = gp_.g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
		statistics[*it] = MMInfo(gp_.g.length(*it) + gp_.g.k());
	}
	stream.reset();
	NewExtendedSequenceMapper<Graph> sm(gp_.g, gp_.index, gp_.kmer_mapper, gp_.g.k() + 1);
	while(!stream.eof()) {
		read_type read;
		stream >> read;
		const Sequence &s_read = read.sequence();
		omnigraph::MappingPath<EdgeId> path = sm.MapSequence(s_read);
		if(path.size() == 1 && path[0].second.initial_range.size() == path[0].second.mapped_range.size()) {
			Range initial_range = path[0].second.initial_range;
			Range mapped_range = path[0].second.mapped_range;
			const Sequence &s_edge = gp_.g.EdgeNucls(path[0].first);
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
}

map<EdgeId, MMInfo> CountStatistics(io::IReader<read_type>& stream) {
	map<EdgeId, MMInfo> statistics;
	CountStatistics(stream, statistics);
	return statistics;
}

map<EdgeId, MMInfo> ParallelCountStatistics(io::ReadStreamVector<io::IReader<read_type>> streams) {
    size_t nthreads = streams.size();
    std::vector<map<EdgeId, MMInfo>> statistics(nthreads);

    #pragma omp parallel for num_threads(nthreads) shared(streams, statistics)
    for (size_t i = 0; i < nthreads; ++i) {
        CountStatistics(streams[i], statistics[i]);
    }

   	for(auto it = statistics[0].begin(); it != statistics[0].end(); ++it) {
   		for(size_t j = 1; j < statistics.size(); j++)
   		for(size_t i = 0; i < statistics[0][it->first].size(); i++) {
   			statistics[0][it->first][i] += statistics[j][it->first][i];
   		}
   	}
   	return statistics[0];
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
	Sequence correct = s_mm.Subseq(0, gp_.g.k()) + Sequence(string(1, nucl)) + s_mm.Subseq(gp_.g.k() + 1, gp_.g.k() * 2 + 1);
	if(!gp_.kmer_mapper.CheckCanRemap(s_mm, correct)) {
		return edge;
	}
	VERIFY(nucl != s_mm[gp_.g.k()]);
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

size_t CorrectAllEdges(map<EdgeId, MMInfo> statistics) {
	size_t res = 0;
	for(auto it = gp_.g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
		if(!gp_.g.RelatedVertices(gp_.g.EdgeStart(*it), gp_.g.EdgeEnd(*it)) && statistics.find(*it) != statistics.end()) {
			res += CorrectEdge(*it, statistics[*it]);
		}
	}
	return res;
}

size_t StopMismatchIteration(io::IReader<read_type>& stream) {
	map<EdgeId, MMInfo> statistics = CountStatistics(stream);
	return CorrectAllEdges(statistics);
}

size_t ParallelStopMismatchIteration(io::ReadStreamVector<io::IReader<read_type>> streams) {
	map<EdgeId, MMInfo> statistics = ParallelCountStatistics(streams);
	return CorrectAllEdges(statistics);
}

public:
	MismatchShallNotPass(graph_pack &gp, double relative_threshold = 2) : gp_(gp), relative_threshold_(relative_threshold) {
		VERIFY(relative_threshold >= 1);
	}


	size_t StopAllMismatches(io::IReader<read_type>& stream, size_t max_iterations = 1) {
		size_t res = 0;
		while(max_iterations > 0) {
			size_t last = StopMismatchIteration(stream);
			res += last;
			if(last == 0)
				break;
			max_iterations--;
		}
		return res;
	}

	size_t ParallelStopAllMismatches(io::ReadStreamVector<io::IReader<read_type>> streams, size_t max_iterations = 1) {
		size_t res = 0;
		while(max_iterations > 0) {
			size_t last = ParallelStopMismatchIteration(streams);
			res += last;
			if(last == 0)
				break;
			max_iterations--;
		}
		return res;
	}
};

}
