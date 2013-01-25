//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef OBSERVABLE_GRAPH_HPP_
#define OBSERVABLE_GRAPH_HPP_

#include "omni_utils.hpp"
#include "id_track_handler.hpp"
//#include <limits>

namespace omnigraph {

//template<typename VertexId, typename EdgeId>
//class ReliableComparator {
//private:
//	const BaseIdTrackHandler<VertexId, EdgeId> *int_ids_;
//
//	template<class Element>
//	int GetFakeIntId(Element a) const {
//		if (a.get() == typename Element::pointer_type(1))
//			return numeric_limits<
//					typename BaseIdTrackHandler<VertexId, EdgeId>::realIdType>::min();
//		if (a.get() == typename Element::pointer_type(-1))
//			return numeric_limits<
//					typename BaseIdTrackHandler<VertexId, EdgeId>::realIdType>::max();
//		return int_ids_->ReturnIntId(a);
//	}
//
//public:
//	ReliableComparator(const BaseIdTrackHandler<VertexId, EdgeId> *int_ids) :
//			int_ids_(int_ids) {
//	}
//
//	bool operator()(VertexId a, VertexId b) const {
//
//		return GetFakeIntId(a) < GetFakeIntId(b);
//	}
//
//	bool operator()(EdgeId a, EdgeId b) const {
//		VERIFY(GetFakeIntId(a) != 0 && GetFakeIntId(b) != 0);
//		return GetFakeIntId(a) < GetFakeIntId(b);
//	}
//
//	template<class Element>
//	bool IsValidId(Element a) const {
//		return int_ids_->ReturnIntId(a) != 0;
//	}
//
//	template<class Element>
//	bool IsAFAKEMin(Element a) const {
//		return a == Element(typename Element::pointer_type(1));
//	}
//
//	template<class Element>
//	bool IsAFAKEMax(Element a) const {
//		return a == Element(typename Element::pointer_type(-1));
//	}
//
//	template<class Element>
//	bool IsAFAKE(Element a) const {
//		return IsAFAKEMin(a) || IsAFAKEMax(a);
//	}
//};

template<typename VertexIdT, typename EdgeIdT, typename VertexIterator/* = typename set<VertexIdT>::iterator*/>
class ObservableGraph: private boost::noncopyable {
public:
	typedef VertexIdT VertexId;
	typedef EdgeIdT EdgeId;
	typedef HandlerApplier<VertexId, EdgeId> Applier;
  typedef typename VertexId::type::edge_const_iterator edge_const_iterator;

//	typedef ReliableComparator<VertexId, EdgeId> Comparator;

//	typedef ReliableComparator<VertexId> ReliableVertexComparator;
//	typedef ReliableComparator<EdgeId> ReliableEdgeComparator;

	typedef SmartVertexIterator<ObservableGraph> SmartVertexIt;
	typedef SmartEdgeIterator<ObservableGraph> SmartEdgeIt;

private:
	typedef ActionHandler<VertexId, EdgeId> Handler;

	const HandlerApplier<VertexId, EdgeId> *applier_;

	mutable vector<Handler*> action_handler_list_;

//public:
//	GraphIdTrackHandler<ObservableGraph> element_order_;

protected:
	virtual void FireAddingVertex(VertexId v) {
		for (auto it = action_handler_list_.begin();
				it != action_handler_list_.end(); ++it) {
			applier_->ApplyAdding(*it, v);
		}
	}

	virtual void FireAddingEdge(EdgeId edge) {
		for (auto it = action_handler_list_.begin();
				it != action_handler_list_.end(); ++it) {
			applier_->ApplyAdding(*it, edge);
		}
	}

	virtual void FireAddVertex(VertexId v) {
		for (auto it = action_handler_list_.begin();
				it != action_handler_list_.end(); ++it) {
			applier_->ApplyAdd(*it, v);
		}
	}

	virtual void FireAddEdge(EdgeId edge) {
		for (auto it = action_handler_list_.begin();
				it != action_handler_list_.end(); ++it) {
			applier_->ApplyAdd(*it, edge);
		}
	}

	virtual void FireDeleteVertex(VertexId v) {
		for (auto it = action_handler_list_.rbegin();
				it != action_handler_list_.rend(); ++it) {
			applier_->ApplyDelete(*it, v);
		}
	}

	virtual void FireDeleteEdge(EdgeId edge) {
		TRACE("FireDeleteEdge for "<<action_handler_list_.size()<<" handlers");
		for (auto it = action_handler_list_.rbegin();
				it != action_handler_list_.rend(); ++it) {
			TRACE("FireDeleteEdge to handler "<<(*it)->name());
			applier_->ApplyDelete(*it, edge);
		}
		TRACE("FireDeleteEdge OK");
	}

	virtual void FireMerge(vector<EdgeId> oldEdges, EdgeId newEdge) {
		TRACE("Fire Merge");
		for (auto it = action_handler_list_.begin();
				it != action_handler_list_.end(); ++it) {
			applier_->ApplyMerge(*it, oldEdges, newEdge);
		}
	}

	virtual void FireGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) {
		TRACE("FireGlue for "<<action_handler_list_.size()<<" handlers");
		for (auto it = action_handler_list_.begin();
				it != action_handler_list_.end(); ++it) {
			TRACE("FireGlue to handler "<<(*it)->name());
			applier_->ApplyGlue(*it, new_edge, edge1, edge2);
		}
		TRACE("FireGlue OK");
	}

	virtual void FireSplit(EdgeId edge, EdgeId newEdge1, EdgeId newEdge2) {
		TRACE("Fire Split");
		for (auto it = action_handler_list_.begin();
				it != action_handler_list_.end(); ++it) {
			applier_->ApplySplit(*it, edge, newEdge1, newEdge2);
		}
	}

public:
	virtual void FireVertexSplit(VertexId newVertex,
			vector<pair<EdgeId, EdgeId> > newEdges,
			vector<double> &split_coefficients, VertexId oldVertex) {
		DEBUG("Fire VertexSplit");
		for (auto it = action_handler_list_.begin();
				it != action_handler_list_.end(); ++it) {
			applier_->ApplyVertexSplit(*it, newVertex, newEdges,
					split_coefficients, oldVertex);
		}
	}

//	void PrintHandlers() {
//		cout << "Printing handlers" << endl;
//		for (auto it = action_handler_list_.begin();
//				it != action_handler_list_.end(); ++it) {
//			cout << (*it)->name() << endl;
//		}
//		cout << "End of handlers" << endl;
//	}

	ObservableGraph(HandlerApplier<VertexId, EdgeId> *applier) :
			applier_(applier)/*, element_order_(*this)*/{
	}

	virtual ~ObservableGraph() {
		TRACE("~ObservableGraph")
		delete applier_;
		TRACE("~ObservableGraph ok")
	}

	void AddActionHandler(Handler* action_handler) const {
		#pragma omp critical(action_handler_list_modification)
		{
			TRACE("Action handler " << action_handler->name() << " added");
			if (find(action_handler_list_.begin(), action_handler_list_.end(),
					action_handler) != action_handler_list_.end()) {
				VERIFY_MSG(false,
						"Action handler " << action_handler->name() << " has already been added");
			} else {
				action_handler_list_.push_back(action_handler);
			}
		}
	}

	bool RemoveActionHandler(const Handler* action_handler) const {
		bool result = false;
		#pragma omp critical(action_handler_list_modification)
		{
			auto it = std::find(action_handler_list_.begin(), action_handler_list_.end(), action_handler);
			if (it != action_handler_list_.end()) {
				action_handler_list_.erase(it);
				TRACE("Action handler " << action_handler->name() << " removed");
				result = true;
			} else {
				TRACE("Action handler " << action_handler->name() << " wasn't found among graph action handlers");
			}
		}

		return result;
	}

	virtual VertexIterator begin() const = 0;

	virtual VertexIterator end() const = 0;

	//todo think of moving to AbstractGraph
	virtual const vector<EdgeId> OutgoingEdges(VertexId vertex) const = 0;

  virtual edge_const_iterator out_begin(VertexId v) const = 0;

  virtual edge_const_iterator out_end(VertexId v) const = 0;

	template<typename Comparator>
	SmartVertexIterator<ObservableGraph, Comparator> SmartVertexBegin(
			const Comparator& comparator) const {
		return SmartVertexIterator<ObservableGraph, Comparator>(*this,
				comparator);
	}

	SmartVertexIterator<ObservableGraph> SmartVertexBegin() const {
		return SmartVertexIterator<ObservableGraph>(*this);
	}

	template<typename Comparator>
	SmartEdgeIterator<ObservableGraph, Comparator> SmartEdgeBegin(
			const Comparator& comparator) const {
		return SmartEdgeIterator<ObservableGraph, Comparator>(*this, comparator);
	}

	SmartEdgeIterator<ObservableGraph> SmartEdgeBegin() const {
		return SmartEdgeIterator<ObservableGraph>(*this);
	}

	//Use very carefully!
	void FireProject(EdgeId edge1, EdgeId edge2) {
		FireGlue(edge2, edge1, edge2);
	}

	const Applier& GetHandlerApplier() const {
		return *applier_;
	}

//	ReliableComparator<VertexId> ReliableComparatorInstance() {
//		return ReliableComparator<VertexId>(element_order_);
//	}

	bool AllHandlersThreadSafe() const {
		BOOST_FOREACH(Handler* handler, action_handler_list_) {
			if (!handler->IsThreadSafe()) {
				return false;
			}
		}
		return true;
	}

	// TODO: for debug. remove.
	void PrintHandlersNames() const {
		BOOST_FOREACH(Handler* handler, action_handler_list_) {
			cout << handler->name() << endl;
		}
	}

private:
	DECL_LOGGER("ObservableGraph")
};
}
#endif /* OBSERVABLE_GRAPH_HPP_ */
