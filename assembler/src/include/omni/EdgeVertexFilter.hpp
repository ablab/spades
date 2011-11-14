/*
 * EdgeFilter.hpp
 *
 *  Created on: 29.07.2011
 *
 */

#ifndef EDGEFILTER_HPP_
#define EDGEFILTER_HPP_

template<class Graph>
class EdgeVertexFilter {
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	set<VertexId> vertices_;
	set<EdgeId> edges_;
	Graph &g_;
public:
	EdgeVertexFilter(Graph &g, const vector<VertexId>& vertices, bool add_verices_with_conjugate): g_(g){
		for (auto it = vertices.begin(); it != vertices.end(); ++it){
			vertices_.insert(*it);
			if (add_verices_with_conjugate){
				vertices_.insert(g_.conjugate(*it));
			}
		}
		for (auto v_it = vertices_.begin(); v_it != vertices_.end(); ++v_it) {
			const vector<EdgeId> edges = g_.OutgoingEdges(*v_it);
			TRACE("working with vertex " << *v_it);
			for (auto e_it = edges.begin(); e_it != edges.end(); ++e_it) {
				VertexId edge_end = g_.EdgeEnd(*e_it);
				TRACE(g_.coverage(*e_it)<<" " << g_.length(*e_it));
				if (vertices_.count(edge_end) > 0) {
					edges_.insert(*e_it);
					TRACE("Edge added");
				}
			}
		}
	}

	size_t VertexCount(){return vertices_.size();};

	bool EdgeIsPresent(EdgeId e_id){
		if (edges_.count(e_id) > 0) return true;
		return false;
	};

	typename set<EdgeId>::iterator EdgesBegin(){
		return edges_.begin();
	}
	typename set<EdgeId>::iterator EdgesEnd(){
		return edges_.end();
	}
	typename set<VertexId>::iterator VerticesBegin(){
		return vertices_.begin();
	}
	typename set<VertexId>::iterator VerticesEnd(){
		return vertices_.end();
	}

};



#endif /* EDGEFILTER_HPP_ */
