/*
 * tip_clipper.hpp
 *
 *  Created on: Mar 25, 2011
 *      Author: sergey
 */

#ifndef TIP_CLIPPER_HPP_
#define TIP_CLIPPER_HPP_

#include <set>

namespace condensed_graph {

using std::set;

class TipClipper {
	CondensedGraph *g_;
	size_t coverage_gap_;

	size_t MaxCoverageOfRightNeighbours(Vertex* v) {
		size_t max = 0;
		for (char i = 0; i < 4; ++i) {
			if (v -> coverage(i) > max) {
				max = v -> coverage(i);
			}
		}
		return max;
	}

	bool IsTip(Vertex* v) {
		if (g->IsFirst(v) && (v->size() < 2 * g->k())
				&& v->RightNeighbourCount() == 1) {
			Vertex* next = *(g->RightNeighbours(v).begin());
			if (v->coverage(next->nucls()[g_->k() - 1]) * coverage_gap
					< MaxCoverageOfRightNeighbours(next -> complement())) {
				return true;
			}
		}
		return false;
	}

public:
	TipClipper(CondensedGraph *g, size_t coverage_gap) :
		g_(g), coverage_gap_(coverage_gap) {

	}
	void ClipTips();
};

void TipClipper::ClipTips() {
	for (set<Vertex*>::const_iterator v_it = g->vertices().begin(); v_it
			!= g->vertices().end(); ++v_it) {
		Vertex* v = *v_it;
		if (IsTip(v)) {
			g_->UnLinkAll(v);
			g_->DeleteVertex(v);
		}
	}
}

}
#endif /* TIP_CLIPPER_HPP_ */
