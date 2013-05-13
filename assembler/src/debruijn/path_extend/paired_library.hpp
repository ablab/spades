//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * paired_library.hpp
 *
 *  Created on: Feb 19, 2012
 *      Author: andrey
 */

#ifndef PAIRED_LIBRARY_HPP_
#define PAIRED_LIBRARY_HPP_

#include "../debruijn_graph.hpp"

//#include "../new_debruijn.hpp"
//#include "../graph_pack.hpp"

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
	bool is_mate_pair_;

	double single_threshold_;
	double coverage_coeff_;

	PairedInfoIndexT<Graph>& index_;
	PairedInfoIndexT<Graph>& index_not_cl_;

	PairedInfoLibrary(size_t k, Graph& g, size_t readS, size_t insS, size_t var,
			PairedIndexT& index,  bool is_mate_pair_ = false) :
			g_(g), k_(k), read_size_(readS), insert_size_(insS), is_variation_(var),
			is_mate_pair_(is_mate_pair_), index_(index), index_not_cl_(index_) {
		setDefault();
	}

	PairedInfoLibrary(size_t k, Graph& g, size_t readS, size_t insS, size_t var,
			PairedIndexT& index, PairedIndexT& index_not_cl, bool is_mate_pair_ = false) :
			g_(g), k_(k), read_size_(readS), insert_size_(insS), is_variation_(var),
			is_mate_pair_(is_mate_pair_), index_(index), index_not_cl_(index_not_cl) {
		setDefault();
	}

	void setDefault(){
		coverage_coeff_ = 1.0;
		single_threshold_ = -1.0;
	}

	PairedInfoLibrary(const PairedInfoLibrary& lib) :
			g_(lib.g_), k_(lib.k_), read_size_(lib.read_size_), insert_size_(lib.insert_size_), is_variation_(lib.is_variation_),
					is_mate_pair_(lib.is_mate_pair_), single_threshold_(lib.single_threshold_), coverage_coeff_(lib.coverage_coeff_),
					index_(lib.index_), index_not_cl_(lib.index_not_cl_) {

	}

    void SetCoverage(double cov) {
        coverage_coeff_ = cov;
    }

    void SetSingleThreshold(double threshold){
    	single_threshold_ = threshold;
    }

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

        int distanceDev = (int) pointIter->var; //max((int) pointIter->var, (int) is_variation_);

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

	double get_all_pi_count(EdgeId e, bool first) {
		vector<PairInfo<EdgeId> > all_pi = index_not_cl_.GetEdgeInfo(e);
		double result = 0;
		for (size_t i = 0; i < all_pi.size(); ++i) {
			if ((all_pi[i].point.d > 0 and first)
					or (all_pi[i].point.d < 0 and !first)) {
				result += all_pi[i].point.weight;
			}
		}
		return result;
	}

	double get_all_norm_pi(EdgeId e, bool first) {
		return get_all_norm_pi_and_size(e, first, index_not_cl_).first;
	}

	double get_all_norm_pi_cl(EdgeId e, bool first) {
		return get_all_norm_pi_and_size(e, first, index_).first;
	}

	pair<double, size_t> get_all_norm_pi_and_size(EdgeId e, bool first, PairedInfoIndexT<Graph>& index) {
			vector<PairInfo<EdgeId> > all_pi = index.GetEdgeInfo(e);
			map<PairInfo<EdgeId>, double> result_pi;
			clust_pi(all_pi, result_pi);
			double result = 0;
			size_t size = 0;
			for (auto iter = result_pi.begin(); iter != result_pi.end(); ++iter) {
				if (math::gr(iter->second, 0.)
						and ((iter->first.point.d > 0 and first)
								or (iter->first.point.d < 0 and !first))) {
					result += iter->first.weight() / iter->second;
					size += 1;
				}
			}
			return make_pair(result, size);
		}

	double get_all_norm_pi_aver(EdgeId e, bool first) {
		pair<double, size_t> res = get_all_norm_pi_and_size(e, first, index_not_cl_);
		return res.first / (double)res.second;
	}

	double get_all_norm_pi_aver_cl(EdgeId e, bool first) {
		pair<double, size_t> res = get_all_norm_pi_and_size(e, first, index_);
		return res.first / (double) res.second;
	}

	void clust_pi(vector<PairInfo<EdgeId> >& pair_edges,
			map<PairInfo<EdgeId>, double>& result) {
		map<pair<EdgeId, EdgeId>, vector<Point> > pi_map;
		for (auto iter = pair_edges.begin(); iter != pair_edges.end(); ++iter) {
			pair<EdgeId, EdgeId> pi(iter->first, iter->second);
			if (pi_map.find(pi) == pi_map.end()) {
				pi_map.insert(make_pair(pi, vector<Point>()));
			}
			pi_map[pi].push_back(iter->point);
		}
		clust_pi(pi_map, result);
	}

	void clust_pi(map<pair<EdgeId, EdgeId>, vector<Point> >& pi_map,
			map<PairInfo<EdgeId>, double>& result) {
		for (auto iter = pi_map.begin(); iter != pi_map.end(); ++iter) {
			std::sort((iter->second).begin(), (iter->second).end());
			double d = iter->second[iter->second.size() / 2].d;
			double w = find_weight(iter->second);
			PairInfo<EdgeId> pi((*iter).first.first, (*iter).first.second,
					Point(d, w, 0));
			result[pi] = IdealPairedInfo((*iter).first.first,
					(*iter).first.second, d);
		}
	}

	double find_weight(vector<Point>& points) {
		double w = 0;
		for (auto point = points.begin(); point != points.end(); ++point) {
			w += (*point).weight;
		}
		return w;
	}
};

typedef std::vector<PairedInfoLibrary *> PairedInfoLibraries;

} // path extend


#endif /* PAIRED_LIBRARY_HPP_ */
