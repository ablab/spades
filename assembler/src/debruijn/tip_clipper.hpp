///*
// * tip_clipper.hpp
// *
// *  Created on: Mar 25, 2011
// *      Author: sergey
// */
//
//#ifndef TIP_CLIPPER_HPP_
//#define TIP_CLIPPER_HPP_
//
//#include <set>
//
//namespace condensed_graph {
//
//using std::set;
//
//class TipClipper {
//public:
//	void ClipTips(CondensedGraph *g);
//};
//
//void TipClipper::ClipTips(CondensedGraph *g) {
//	for (set<Vertex*>::const_iterator v_it = g->vertices().begin(); v_it
//			!= g->vertices().begin(); ++v_it) {
//		Vertex* v = *v_it;
//		if (v->size() < 2*g->k() && v->RightNeighbourCount() == 1) {
//			*(g->RightNeighbours(v).begin())
//		}
//	}
//}
//
//}
//#endif /* TIP_CLIPPER_HPP_ */
