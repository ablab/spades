//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "openmp_wrapper.h"
#include "standard.hpp"

#include "io/paired_read.hpp"
#include "omni/omni_utils.hpp"
#include "omni/visualization/graph_colorer.hpp"
#include "omni/id_track_handler.hpp"
#include "omni/splitters.hpp"
#include "omni/path_processor.hpp"

#include "logger/logger.hpp"
#include "xmath.h"
#include "sequence/sequence_tools.hpp"

#include "runtime_k.hpp"

#include "path_helper.hpp"

#include "debruijn_graph.hpp"
#include "indices/perfect_hash_map.hpp"
#include "edge_index.hpp"

#include <iostream>

namespace debruijn_graph {

using omnigraph::Path;
using omnigraph::MappingPath;
using omnigraph::Range;
using omnigraph::MappingRange;

inline double PairedReadCountWeight(const MappingRange&, const MappingRange&) {
	return 1.;
}

inline double KmerCountProductWeight(const MappingRange& mr1,
                                     const MappingRange& mr2) {
    return (double)(mr1.initial_range.size() * mr2.initial_range.size());
}

class WeightDEWrapper {
private:

    vector<double> new_hist;
    int left_x;
    int insert_size;

	void ExtendLinear(const std::map<int, size_t> & hist) {
        size_t sum_weight = 0;

        for (auto iter = hist.begin(); iter != hist.end(); ++iter)
            sum_weight += iter->second;
        DEBUG(sum_weight);

        VERIFY(hist.size() > 0);
        auto iter = hist.begin();

        left_x = iter->first;

        int prev = iter->first;
        size_t prev_val = iter->second;

        new_hist.push_back((double)prev_val / (double)sum_weight);
        ++iter;

        for (; iter != hist.end(); ++iter) {
            int x = iter->first;
            size_t y = iter->second;
            double tan = ((double)y - (double)prev_val) / (x - prev);

            VERIFY(prev < x);
            for (int i = prev + 1; i <= x; ++i) {
                new_hist.push_back(((double)prev_val + tan * (i - prev)) / (double)sum_weight);
            }
            prev = x;
            prev_val = y;
            DEBUG("hist " << x << " " << y);
        }
	}

public:
    WeightDEWrapper(const map<int, size_t>& hist, double IS) {
        DEBUG("WeightDEWrapper " << IS);
        insert_size = (int) IS;
        DEBUG("Extending linear");
        ExtendLinear(hist);
    }

    ~WeightDEWrapper() {
    }


    double CountWeight(int x) const {
        int xx = insert_size - left_x + x - 1;

        if (!(xx >= 0 && xx < (int) new_hist.size())) return 0.;
        VERIFY(math::le(new_hist[xx], 1.));
        return 1000. * new_hist[xx];
    }
};

inline double UnityFunction(int /*x*/) {
    return 1.;
}

template<class Graph>
Sequence MergeSequences(const Graph& g,
		const vector<typename Graph::EdgeId>& continuous_path) {
	vector < Sequence > path_sequences;
	path_sequences.push_back(g.EdgeNucls(continuous_path[0]));
	for (size_t i = 1; i < continuous_path.size(); ++i) {
		VERIFY(
				g.EdgeEnd(continuous_path[i - 1])
						== g.EdgeStart(continuous_path[i]));
		path_sequences.push_back(g.EdgeNucls(continuous_path[i]));
	}
	return MergeOverlappingSequences(path_sequences, g.k());
}

}
