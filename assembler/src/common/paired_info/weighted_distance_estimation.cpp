#include "weighted_distance_estimation.hpp"

namespace omnigraph {
namespace de {

using namespace debruijn_graph;

WeightedDistanceEstimator::EstimHist WeightedDistanceEstimator::EstimateEdgePairDistances(EdgePair ep,
                                                                                          const InHistogram &histogram,
                                                                                          const GraphLengths &raw_forward) const {
    using std::abs;
    using namespace math;
    TRACE("Estimating with weight function");
    size_t first_len = this->graph().length(ep.first);
    size_t second_len = this->graph().length(ep.second);

    EstimHist result;
    int maxD = rounded_d(histogram.max()), minD = rounded_d(histogram.min());
    vector<int> forward;
    for (auto len : raw_forward) {
        int length = (int) len;
        if (minD - (int) this->max_distance_ <= length && length <= maxD + (int) this->max_distance_) {
            forward.push_back(length);
        }
    }
    if (forward.size() == 0)
        return result;

    DEDistance max_dist = this->max_distance_;
    size_t i = 0;
    vector<double> weights(forward.size());
    for (auto point : histogram) {
        DEDistance cur_dist(forward[i]), next_dist(forward[i + 1]);
        if (le(2 * point.d + DEDistance(second_len), DEDistance(first_len)))
            continue;
        while (i + 1 < forward.size() && next_dist < point.d) {
            ++i;
        }
        if (i + 1 < forward.size() && ls(DEDistance(next_dist) - point.d, point.d - DEDistance(cur_dist))) {
            ++i;
            if (le(abs(cur_dist - point.d), max_dist))
                weights[i] += point.weight * weight_f_(forward[i] - rounded_d(point));
        }
        else if (i + 1 < forward.size() && eq(next_dist - point.d, point.d - cur_dist)) {
            if (le(abs(cur_dist - point.d), max_dist))
                weights[i] += point.weight * 0.5 * weight_f_(forward[i] - rounded_d(point));

            ++i;

            if (le(abs(cur_dist - point.d), max_dist))
                weights[i] += point.weight * 0.5 * weight_f_(forward[i] - rounded_d(point));
        } else if (le(abs(cur_dist - point.d), max_dist))
            weights[i] += point.weight * weight_f_(forward[i] - rounded_d(point));
    }

    for (size_t i = 0; i < forward.size(); ++i)
        if (gr(weights[i], 0.))
            result.push_back(make_pair(forward[i], weights[i]));

    return result;
}
}
}
