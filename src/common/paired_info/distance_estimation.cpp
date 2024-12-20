//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "distance_estimation.hpp"
#include "pair_info_bounds.hpp"
#include "assembly_graph/paths/path_processor.hpp"

namespace omnigraph::de {

using namespace debruijn_graph;

std::vector<size_t> GraphDistanceFinder::GetGraphDistancesLengths(EdgeId e1, EdgeId e2) const {
    LengthMap m;
    m.insert({e2, {}});

    FillGraphDistancesLengths(e1, m);

    return m[e2];
}

void GraphDistanceFinder::FillGraphDistancesLengths(EdgeId e1, LengthMap &second_edges) const {
    std::vector<size_t> path_lower_bounds;
    size_t path_upper_bound = PairInfoPathLengthUpperBound(graph_.k(), insert_size_, delta_);
    PathProcessor <Graph> paths_proc(graph_, graph_.EdgeEnd(e1), path_upper_bound);

    for (auto &entry : second_edges) {
        EdgeId e2 = entry.first;
        size_t path_lower_bound = PairInfoPathLengthLowerBound(graph_.k(), graph_.length(e1),
                                                               graph_.length(e2), gap_, delta_);

        TRACE("Bounds for paths are " << path_lower_bound << " " << path_upper_bound);

        DistancesLengthsCallback<Graph> callback(graph_);
        paths_proc.Process(graph_.EdgeStart(e2), path_lower_bound, path_upper_bound, callback);
        GraphLengths lengths = callback.distances();
        for (size_t j = 0; j < lengths.size(); ++j) {
            lengths[j] += graph_.length(e1);
            TRACE("Resulting distance set for " <<
                                                " edge " << graph_.int_id(e2) <<
                                                " #" << j << " length " << lengths[j]);
        }

        if (e1 == e2)
            lengths.push_back(0);

        std::sort(lengths.begin(), lengths.end());
        entry.second = lengths;
    }
}

void AbstractDistanceEstimator::FillGraphDistancesLengths(EdgeId e1, LengthMap &second_edges) const {
    distance_finder_.FillGraphDistancesLengths(e1, second_edges);
}

AbstractDistanceEstimator::OutHistogram AbstractDistanceEstimator::ClusterResult(EdgePair,
                                                                                 const EstimHist &estimated) const {
    OutHistogram result;
    for (size_t i = 0; i < estimated.size(); ++i) {
        size_t left = i;
        DEWeight weight = DEWeight(estimated[i].second);
        while (i + 1 < estimated.size() &&
               (estimated[i + 1].first - estimated[i].first) <= (int) linkage_distance_) {
            ++i;
            weight += estimated[i].second;
        }
        DEDistance center = DEDistance((estimated[left].first + estimated[i].first) * 0.5);
        DEVariance var = DEVariance((estimated[i].first - estimated[left].first) * 0.5);
        result.insert(Point(center, weight, var));
    }
    return result;
}

void DistanceEstimator::Estimate(PairedInfoIndexT<Graph> &result, size_t nthreads) const  {
    this->Init();
    const auto &index = this->index();
    ConcurrentUnorderedClusteredPairedInfoBuffer<Graph> buffer(graph());

    omnigraph::IterationHelper<Graph, EdgeId> edges(graph());
    auto ranges = edges.Ranges(nthreads * 16);

    DEBUG("Processing");
#   pragma omp parallel for schedule(guided) num_threads(nthreads)
    for (size_t i = 0; i < ranges.size(); ++i) {
        TRACE("Processing chunk #" << i);

        for (EdgeId e : ranges[i]) {
            TRACE("Estimating for edge " << e);
            ProcessEdge(e, index, buffer);
        }
    }

    result.Merge(buffer);
}

DistanceEstimator::EstimHist DistanceEstimator::EstimateEdgePairDistances(EdgePair ep, const InHistogram &histogram,
                                                                          const GraphLengths &raw_forward) const {
    using std::abs;
    using namespace math;
    EdgeId e1 = ep.first, e2 = ep.second;
    size_t first_len = this->graph().length(e1), second_len = this->graph().length(e2);
    int minD = rounded_d(histogram.min()), maxD = rounded_d(histogram.max());

    TRACE("Bounds are " << minD << " " << maxD);
    EstimHist result;
    std::vector<DEDistance> forward;
    forward.reserve(raw_forward.size());
    for (auto raw_length : raw_forward) {
        int length = int(raw_length);
        if (minD - int(max_distance_) <= length && length <= maxD + int(max_distance_))
            forward.push_back(DEDistance(length));
    }
    if (forward.size() == 0)
        return result;

    size_t cur_dist = 0;
    std::vector<DEWeight> weights(forward.size(), 0);
    for (auto point : histogram) {
        if (ls(2 * point.d + DEDistance(second_len), DEDistance(first_len)))
            continue;
        while (cur_dist + 1 < forward.size() && forward[cur_dist + 1] < point.d)
            ++cur_dist;

        if (cur_dist + 1 < forward.size() &&
            ls(forward[cur_dist + 1] - point.d, point.d - forward[cur_dist])) {
            ++cur_dist;

            if (le(abs(forward[cur_dist] - point.d), max_distance_))
                weights[cur_dist] += point.weight;
        } else if (cur_dist + 1 < forward.size() &&
                   eq(forward[cur_dist + 1] - point.d, point.d - forward[cur_dist])) {
            if (le(abs(forward[cur_dist] - point.d), max_distance_))
                weights[cur_dist] += point.weight * 0.5;
            ++cur_dist;
            if (le(abs(forward[cur_dist] - point.d), max_distance_))
                weights[cur_dist] += point.weight * 0.5;
        } else {
            if (le(abs(forward[cur_dist] - point.d), max_distance_))
                weights[cur_dist] += point.weight;
        }
    }

    for (size_t i = 0; i < forward.size(); ++i)
        if (ge(weights[i], DEWeight(0)))
            result.emplace_back(forward[i], weights[i]);

    VERIFY(result.size() == forward.size());
    return result;
}

void DistanceEstimator::ProcessEdge(EdgeId e1, const InPairedIndex &pi, Buffer &result) const {
    typename base::LengthMap second_edges;
    auto inner_map = pi.GetHalf(e1);
    for (auto i : inner_map)
        second_edges[i.first];

    this->FillGraphDistancesLengths(e1, second_edges);

    for (const auto &entry: second_edges) {
        EdgeId e2 = entry.first;
        EdgePair ep(e1, e2);

        VERIFY(ep <= pi.ConjugatePair(ep));

        const GraphLengths &forward = entry.second;
        TRACE("Edge pair is " << this->graph().int_id(ep.first)
                              << " " << this->graph().int_id(ep.second));
        auto hist = pi.Get(e1, e2);
        const EstimHist &estimated = this->EstimateEdgePairDistances(ep, hist, forward);
        OutHistogram res = this->ClusterResult(ep, estimated);
        this->AddToResult(res, ep, result);
    }
}

} // namespace omnigraph::de
