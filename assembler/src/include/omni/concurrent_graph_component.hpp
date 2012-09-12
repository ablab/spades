//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************


/*
 * graph_component_wrapper.hpp
 *
 *  Created on: Aug 10, 2012
 *      Author: Alexander Opeykin (alexander.opeykin@gmail.com)
 */


#ifndef CONCURRENT_GRAPH_COMPONENT_HPP_
#define CONCURRENT_GRAPH_COMPONENT_HPP_


#include <boost/foreach.hpp>

#include "standard_base.hpp"
#include "order_and_law.hpp"
#include "abstract_editable_graph.hpp"


namespace omnigraph {


template <typename Graph>
class ConcurrentGraphComponent
		: public AbstractEditableGraph<
		  	  typename Graph::VertexId,
		  	  typename Graph::EdgeId,
		  	  typename Graph::DataMaster,
		  	  typename unordered_set<typename Graph::VertexId>::const_iterator > {

public:

	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::EdgeData EdgeData;
	typedef typename Graph::VertexData VertexData;
	typedef typename unordered_set<VertexId>::const_iterator VertexIterator;
	typedef typename Graph::DataMaster DataMaster;
	typedef AbstractEditableGraph<VertexId, EdgeId, DataMaster, VertexIterator> base;


	template<class InputVertexIterator>
	ConcurrentGraphComponent(
			Graph& graph,
			HandlerApplier<VertexId, EdgeId>* applier,
			restricted::IdDistributor& id_distributor,
			InputVertexIterator verticesBegin,
			InputVertexIterator verticesEnd)
				: base(applier, graph.master()),
				graph_(graph), vertices_(verticesBegin, verticesEnd),
				edge_id_distributor_(id_distributor) {

		BOOST_FOREACH(const VertexId& vertex, vertices_) {
			if (!IsInComponent(graph_.OutgoingEdges(vertex)) ||
						!IsInComponent(graph_.IncomingEdges(vertex))) {

				border_vertices_.insert(vertex);
			}
		}

		all_actions_valid_ = true;
	}


	virtual const EdgeData& data(EdgeId edge) const {
		SetFlagIfNotInComponent(edge);
		return graph_.data(edge);
	}

	virtual const VertexData& data(VertexId vertex) const {
		SetFlagIfNotInComponent(vertex);
		return graph_.data(vertex);
	}

	virtual size_t OutgoingEdgeCount(VertexId vertex) const {
		return OutgoingEdges(vertex).size();
	}

	virtual size_t IncomingEdgeCount(VertexId vertex) const {
		return IncomingEdges(vertex).size();
	}

	virtual const vector<EdgeId> OutgoingEdges(VertexId vertex) const {
		if (IsInComponent(vertex)) {
			return GetEdgesFromComponentOrSetFlag(graph_.OutgoingEdges(vertex));
		} else {
			// edges can be filtered by GetEdgesFromComponentOrSetFlag.
			// but this execution branch prevents data races.
			TRACE("Invalidate component action on OutgoingEdges for " << str(vertex));
			all_actions_valid_ = false;
			return vector<EdgeId>();
		}
	}

	virtual const vector<EdgeId> IncomingEdges(VertexId vertex) const {
		if (IsInComponent(vertex)) {
			return GetEdgesFromComponentOrSetFlag(graph_.IncomingEdges(vertex));
		} else {
			TRACE("Invalidate component action on IncomingEdges for " << str(vertex));
			all_actions_valid_ = false;
			return vector<EdgeId>();
		}
	}

	virtual vector<EdgeId> GetEdgesBetween(VertexId vertex1, VertexId vertex2) const {
		if (IsInComponent(vertex1) && IsInComponent(vertex2)) {
			return graph_.GetEdgesBetween(vertex1, vertex2);
		} else {
			TRACE("Invalidate component action on GetEdgesBetween for "
					<< str(vertex1) << " and " << str(vertex2));
			all_actions_valid_ = false;
			return vector<EdgeId>();
		}
	}

	virtual VertexId EdgeStart(EdgeId edge) const {
		SetFlagIfNotInComponent(edge);
		return graph_.EdgeStart(edge);
	}

	virtual VertexId EdgeEnd(EdgeId edge) const {
		SetFlagIfNotInComponent(edge);
		return graph_.EdgeEnd(edge);
	}

	virtual bool RelatedVertices(VertexId vertex1, VertexId vertex2) const {
		SetFlagIfNotInComponent(vertex1);
		SetFlagIfNotInComponent(vertex2);
		return graph_.RelatedVertices(vertex1, vertex2);
	}

	virtual bool CanCompressVertex(const VertexId& vertex) const {
		SetFlagIfNotInternal(vertex);
		return graph_.CanCompressVertex(vertex);
	}

	size_t k() const {
		return graph_.k();
	}

	double coverage(EdgeId edge) const {
		SetFlagIfNotInComponent(edge);
		return graph_.coverage(edge);
	}

	size_t length(EdgeId edge) const {
		SetFlagIfNotInComponent(edge);
		return graph_.length(edge);
	}

	size_t length(VertexId vertex) const {
		SetFlagIfNotInComponent(vertex);
		return graph_.length(vertex);
	}

	const Sequence& EdgeNucls(EdgeId edge) const {
		SetFlagIfNotInComponent(edge);
		return graph_.EdgeNucls(edge);
	}


	virtual std::string str(const EdgeId edge) const {
		return graph_.str(edge);
	}

	virtual std::string str(const VertexId vertex) const {
			return graph_.str(vertex);
	}

	template<typename Comparator>
	SmartEdgeIterator<ConcurrentGraphComponent, Comparator> SmartEdgeBegin(
			const Comparator& comparator, vector<EdgeId>* edges = 0) const {

		bool all_actions_was_valid = all_actions_valid_;

		auto it = SmartEdgeIterator<ConcurrentGraphComponent, Comparator>(*this, comparator, edges);

		// No doubt that creating of iterator will invalidate component.
		// But it's ok.
		all_actions_valid_ = all_actions_was_valid;

		return it;
	}

	SmartEdgeIterator<ConcurrentGraphComponent> SmartEdgeBegin(vector<EdgeId>* edges = 0) const {
		bool all_actions_was_valid = all_actions_valid_;

		auto it = SmartEdgeIterator<ConcurrentGraphComponent, std::less<EdgeId>>(
				*this, std::less<EdgeId>(), edges);

		// No doubt that creating of iterator will invalidate component.
		// But it's ok.
		all_actions_valid_ = all_actions_was_valid;

		return it;
	}

	virtual VertexIterator begin() const {
		return vertices_.begin();
	}

	virtual VertexIterator end() const {
		return vertices_.end();
	}

	virtual ~ConcurrentGraphComponent() {
		// failing here means that algorithm performed on this component
		// created not temporary vertex (did not delete it)
		VERIFY(temporary_vertices_.size() == 0);
	}

	void ValidateComponent() {
		all_actions_valid_ = true;
	}

	bool IsValid() const {
		return all_actions_valid_;
	}

	void GetEdgesGoingOutOfComponent(vector<EdgeId>& output) {
		BOOST_FOREACH(const VertexId& vertex, vertices_) {
			vector<EdgeId> edges = graph_.OutgoingEdges(vertex);
			BOOST_FOREACH(const EdgeId& edge, edges) {
				if (!IsInComponent(edge)) {
					output.push_back(edge);
				}
			}
		}
	}


// Self methods
	bool IsInComponent(const VertexId& vertex) const {
		return vertices_.find(vertex) != vertices_.end();
	}

	bool IsInComponent(const EdgeId& edge) const {
		return IsInComponent(graph_.EdgeStart(edge)) && IsInComponent(graph_.EdgeEnd(edge));
	}

	bool IsInComponent(const std::vector<VertexId>& vertices) const {
		BOOST_FOREACH(const VertexId& vertex, vertices) {
			if (!IsInComponent(vertex)) {
				return false;
			}
		}

		return true;
	}

	bool IsInComponent(const std::vector<EdgeId>& edges) const {
		BOOST_FOREACH(const EdgeId& edge, edges) {
			if (!IsInComponent(edge)) {
				return false;
			}
		}

		return true;
	}

	bool IsAtBorder(const VertexId& vertex) const {
		return border_vertices_.find(vertex) != border_vertices_.end();
	}

	bool IsInternal(const VertexId& vertex) const {
		return IsInComponent(vertex) && !IsAtBorder(vertex);
	}

	virtual bool IsInternalSafe(const VertexId& vertex) const = 0;

	virtual bool IsInternalSafe(const EdgeId& edge) const = 0;


protected:

	virtual void AddVertexToComponent(VertexId vertex) = 0;

	virtual void DeleteVertexFromComponent(VertexId vertex) = 0;

	virtual bool AdditionalCompressCondition(VertexId vertex) const  {
		SetFlagIfNotInComponent(vertex);
		return graph_.AdditionalCompressCondition(vertex);
	}

	virtual EdgeId HiddenAddEdge(VertexId vertex1, VertexId vertex2, const EdgeData &data) {
		return HiddenAddEdge(vertex1, vertex2, data, &edge_id_distributor_);
	}

	virtual EdgeId HiddenAddEdge(VertexId vertex1, VertexId vertex2,
			const EdgeData &data, restricted::IdDistributor * id_distributor) {
		VERIFY(IsInComponent(vertex1));
		VERIFY(IsInComponent(vertex2));
		VERIFY(all_actions_valid_);
		return graph_.HiddenAddEdge(vertex1, vertex2, data, id_distributor);
	}

	virtual void HiddenDeleteEdge(EdgeId edge) {
		VERIFY(IsInComponent(edge));
		VERIFY(all_actions_valid_);
		graph_.HiddenDeleteEdge(edge);
	}

	virtual vector<EdgeId> CorrectMergePath(const vector<EdgeId>& path) {
		VERIFY(IsInComponent(path)); // TODO: debug only??
		vector<EdgeId> corrected_path = graph_.CorrectMergePath(path);
		VERIFY(IsInComponent(corrected_path)); // TODO: are always from same component??
		return corrected_path;
	}

	virtual vector<EdgeId> EdgesToDelete(const vector<EdgeId> &path) {
		vector<EdgeId> edges_to_delete = graph_.EdgesToDelete(path);
		VERIFY(IsInComponent(edges_to_delete));
		return edges_to_delete;
	}

	virtual vector<VertexId> VerticesToDelete(const vector<EdgeId> &path) {
		vector<VertexId> vertices_to_delete = graph_.VerticesToDelete(path);
		VERIFY(IsInComponent(vertices_to_delete));
		return vertices_to_delete;
	}

	virtual VertexId CreateVertex(const VertexData &data) {
		return graph_.CreateVertex(data);
	}

	virtual void DestroyVertex(VertexId vertex) {
		graph_.DestroyVertex(vertex);
	}

protected:
	// observable graph methods.
	virtual void FireAddingVertex(VertexId vertex) {
		base::FireAddingVertex(vertex);
		graph_.FireAddingVertex(vertex);
	}

	virtual void FireAddingEdge(EdgeId edge) {
		base::FireAddingEdge(edge);
		graph_.FireAddingEdge(edge);
	}

	virtual void FireAddVertex(VertexId vertex) {
		base::FireAddVertex(vertex);
		graph_.FireAddVertex(vertex);
	}

	virtual void FireAddEdge(EdgeId edge) {
		base::FireAddEdge(edge);
		graph_.FireAddEdge(edge);
	}

	virtual void FireDeleteVertex(VertexId vertex) {
		base::FireDeleteVertex(vertex);
		graph_.FireDeleteVertex(vertex);
	}

	virtual void FireDeleteEdge(EdgeId edge) {
		base::FireDeleteEdge(edge);
		graph_.FireDeleteEdge(edge);
	}

	virtual void FireMerge(vector<EdgeId> oldEdges, EdgeId newEdge) {
		base::FireMerge(oldEdges, newEdge);
		graph_.FireMerge(oldEdges, newEdge);
	}

	virtual void FireGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) {
		base::FireGlue(new_edge, edge1, edge2);
		graph_.FireGlue(new_edge, edge1, edge2);
	}

	virtual void FireVertexSplit(VertexId newVertex,
				vector<pair<EdgeId, EdgeId> > newEdges,
				vector<double> &split_coefficients, VertexId oldVertex) {
		base::FireVertexSplit(newVertex, newEdges, split_coefficients, oldVertex);
		graph_.FireVertexSplit(newVertex, newEdges, split_coefficients, oldVertex);
	}


	const vector<EdgeId> GetEdgesFromComponentOrSetFlag(const std::vector<EdgeId>& edges) const {
		vector<EdgeId> edges_from_component;
		edges_from_component.reserve(edges.size());

		BOOST_FOREACH(const EdgeId& edge, edges) {
			if (IsInComponent(edge)) {
				edges_from_component.push_back(edge);
			} else {
				all_actions_valid_ = false;
				TRACE("COMPONENT INVALIDATION: Edge Not in component " << graph_.str(edge));
			}
		}

		return edges_from_component;
	}

	void SetFlagIfNotInternal(const VertexId& vertex) const {
		if (!IsInternal(vertex)) {
			all_actions_valid_ = false;
			TRACE("COMPONENT INVALIDATION: vertex not internal " << graph_.str(vertex));
		}
	}

	void SetFlagIfNotInComponent(const EdgeId& edge) const {
		if (!IsInComponent(edge)) {
			all_actions_valid_ = false;
			TRACE("COMPONENT INVALIDATION: Edge Not in component " << graph_.str(edge));
		}
	}

	void SetFlagIfNotInComponent(const VertexId& vertex) const {
		if (!IsInComponent(vertex)) {
			all_actions_valid_ = false;
			TRACE("COMPONENT INVALIDATION: vertex not internal " << graph_.str(vertex));
		}
	}


protected:

	Graph& graph_;

	unordered_set<VertexId> vertices_;
	unordered_set<VertexId> border_vertices_;

	// used for temporary vertices created by algorithms
	unordered_set<VertexId> temporary_vertices_;

	restricted::IdDistributor& edge_id_distributor_;
	mutable bool all_actions_valid_;


private:
	DECL_LOGGER("ConcurrentGraphComponent");
};


} //namespace omnigraph



#endif /* CONCURRENT_GRAPH_COMPONENT_HPP_ */
