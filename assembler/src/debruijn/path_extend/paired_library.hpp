/*
 * paired_library.hpp
 *
 *  Created on: Feb 19, 2012
 *      Author: andrey
 */

#ifndef PAIRED_LIBRARY_HPP_
#define PAIRED_LIBRARY_HPP_

#include "../new_debruijn.hpp"
#include "xmath.h"

using debruijn_graph::Graph;
using debruijn_graph::EdgeId;
using debruijn_graph::VertexId;

using omnigraph::PairedInfoIndex;
using omnigraph::PairInfo;

namespace path_extend {

struct PairedInfoLibrary {

    Graph& g_;

    size_t k_;
    size_t read_size_;
    size_t insert_size_;
    size_t is_variation_;

    double coverage_coeff_;

    PairedInfoIndexT<Graph>& index_;

    PairedInfoLibrary(size_t k, Graph& g, size_t readS, size_t insS, size_t var, PairedInfoIndexT<Graph>& index):
        g_(g), k_(k), read_size_(readS), insert_size_(insS), is_variation_(var), index_(index) {

        coverage_coeff_ = 1.0;
    }

    PairedInfoLibrary(const PairedInfoLibrary& lib): g_(lib.g_), read_size_(lib.read_size_), insert_size_(lib.insert_size_), is_variation_(lib.is_variation_), index_(lib.index_) {

    }

    void SetCoverage(double cov) {
        coverage_coeff_ = cov;
    }


//    void GapStat() {
//        for (auto iter = g_.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
//            int total_sources = 0;
//            int sources_wpi = 0;
//            int sources_wpi_to_sinks = 0;
//            int sources_wuniquepi_to_sinks = 0;
//
//            if (g_.OutgoingEdgeCount(g_.EdgeEnd(*iter)) == 0) {
//                ++total_sources;
//
//                auto pi = index_.GetEdgeInfo(*iter);
//                for (size_t i = 0; i < pi.size(); ++i) {
//                    if (pi[i].first == *iter && pi[i].second != *iter && pi[i].d > 0) {
//
//                    }
//                }
//
//                if (pi.size() > 0)
//                    ++sources_wpi;
//
//
//            }
//
//            INFO("Graph has " << total_sources << ", of which " << sources_wpi  << " have paired info")
//            INFO(sources_wpi_to_sinks << " have paired info to source, of which " << sources_wuniquepi_to_sinks << " are unique");
//        }
//    }

    set<EdgeId> GetEdges(EdgeId e) {
        set<EdgeId> res;
        typedef set<Point> Histogram;
        const InnerMap<Graph>& pairs = index_.GetEdgeInfo(e, 0); // map[second_edge -> histogram]
        for (auto pairIter = pairs.begin(); pairIter != pairs.end(); ++pairIter)
          res.insert(pairIter->first);
        return res;
    }

    void CountDistances(EdgeId e1, EdgeId e2, vector<int>& dist, vector<double>& w) {
      typedef set<Point> Histogram;
    	if (e1 != e2) {
        Histogram histogram = index_.GetEdgePairInfo(e1, e2);
        for (auto pointIter = histogram.begin(); pointIter != histogram.end(); ++pointIter) {
          int pairedDistance = rounded_d(*pointIter);
          if (pairedDistance >= 0) {
            dist.push_back(pairedDistance);
            w.push_back(pointIter->weight);
          }
        }
    	}
	}

    double CountPairedInfo(EdgeId e1, EdgeId e2, int distance) {
      typedef set<Point> Histogram;
      double weight = 0.0;
      Histogram pairs =  index_.GetEdgePairInfo(e1, e2);

      for (auto pointIter = pairs.begin(); pointIter != pairs.end(); ++pointIter) {
        int pairedDistance = rounded_d(*pointIter);

        int distanceDev = max((int) pointIter->var, (int) is_variation_);

        //Can be modified according to distance comparison
        if (pairedDistance >= distance - distanceDev  && pairedDistance <= distance + distanceDev) {
          weight += pointIter->weight;
        }

      }
      return weight;
    }


    double IdealPairedInfo(EdgeId e1, EdgeId e2, int distance) {
        double w = 0.;
        if (distance == 0 && e1 == e2) {
            w = 0. + g_.length(e1) - insert_size_ + 2 * read_size_ + 1 - k_;
        }
        else {
            if (distance < 0) {
                EdgeId tmp = e1;
                e1 = e2;
                e2 = tmp;
                distance = -distance;
            }
            int gap_len = distance - g_.length(e1);
            int right = std::min(insert_size_, gap_len + g_.length(e2) + read_size_);
            int left = std::max(gap_len, int(insert_size_) - int(read_size_) - int(g_.length(e1)));
            w = 0. + right - left + 1 - k_ + is_variation_;
        }
        return math::gr(w, 0.0) ? w : 0.0;
    }

    double NormalizeWeight(const PairInfo<EdgeId>& pair_info) {
        double w = IdealPairedInfo(pair_info.first, pair_info.second, rounded_d(pair_info));

        double result_weight = pair_info.weight();
        if (math::gr(w, 0.)) {
            result_weight /= w;
        } else {
            result_weight = 0.0;
        }

        return result_weight;
    }

    double normalizeByCoverage(double weight) {
        VERIFY(coverage_coeff_ != 0);
        return weight / coverage_coeff_;
    }
};

typedef std::vector<PairedInfoLibrary *> PairedInfoLibraries;

} // path extend


#endif /* PAIRED_LIBRARY_HPP_ */
