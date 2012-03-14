#ifndef OBSERVABLE_GRAPH_HPP_
#define OBSERVABLE_GRAPH_HPP_

#include "omni_utils.hpp"
#include "id_track_handler.hpp"
//#include <limits>

namespace omnigraph {

template<typename VertexId, typename EdgeId>
class ReliableComparator {
private:
	const BaseIdTrackHandler<VertexId, EdgeId> *int_ids_;

	template<class Element>
	int GetFakeIntId(Element a) const {
		if (a.get() == typename Element::pointer_type(1))
			return numeric_limits<
					typename BaseIdTrackHandler<VertexId, EdgeId>::realIdType>::min();
		if (a.get() == typename Element::pointer_type(-1))
			return numeric_limits<
					typename BaseIdTrackHandler<VertexId, EdgeId>::realIdType>::max();
		return int_ids_->ReturnIntId(a);
	}

public:
	ReliableComparator(const BaseIdTrackHandler<VertexId, EdgeId> *int_ids) :
			int_ids_(int_ids) {
	}

	bool operator()(VertexId a, VertexId b) const {

		return GetFakeIntId(a) < GetFakeIntId(b);
	}

	bool operator()(EdgeId a, EdgeId b) const {
		VERIFY(GetFakeIntId(a) != 0 && GetFakeIntId(b) != 0);
		return GetFakeIntId(a) < GetFakeIntId(b);
	}

	template<class Element>
	bool IsValidId(Element a) const {
		return int_ids_->ReturnIntId(a) != 0;
	}

	template<class Element>
	bool IsAFAKEMin(Element a) const {
		return a == Element(typename Element::pointer_type(1));
	}

	template<class Element>
	bool IsAFAKEMax(Element a) const {
		return a == Element(typename Element::pointer_type(-1));
	}

	template<class Element>
	bool IsAFAKE(Element a) const {
		return IsAFAKEMin(a) || IsAFAKEMax(a);
	}
};

template<typename VertexIdT, typename EdgeIdT, typename VertexIterator/* = typename set<VertexIdT>::iterator*/>
class ObservableGraph : private boost::noncopyable {
public:
	typedef VertexIdT VertexId;
	typedef EdgeIdT EdgeId;

	typedef ReliableComparator<VertexId, EdgeId> Comparator;

//	typedef ReliableComparator<VertexId> ReliableVertexComparator;
//	typedef ReliableComparator<EdgeId> ReliableEdgeComparator;

	typedef SmartVertexIterator<ObservableGraph, ObservableGraph::Comparator> SmartVertexIt;
	typedef SmartEdgeIterator<ObservableGraph, ObservableGraph::Comparator> SmartEdgeIt;

private:
	typedef ActionHandler<VertexId, EdgeId> Handler;

	const HandlerApplier<VertexId, EdgeId> *applier_;

	mutable vector<Handler*> action_handler_list_;

public:
	GraphIdTrackHandler<ObservableGraph> element_order_;

protected:
	void FireAddingVertex(VertexId v) {
		for (auto it = action_handler_list_.begin(); it
				!= action_handler_list_.end(); ++it) {
			applier_->ApplyAdding(*it, v);
		}
	}

	void FireAddingEdge(EdgeId edge) {
		for (auto it = action_handler_list_.begin(); it
				!= action_handler_list_.end(); ++it) {
			applier_->ApplyAdding(*it, edge);
		}
	}

	void FireAddVertex(VertexId v) {
		for (auto it = action_handler_list_.begin(); it
				!= action_handler_list_.end(); ++it) {
			applier_->ApplyAdd(*it, v);
		}
	}

	void FireAddEdge(EdgeId edge) {
		for (auto it = action_handler_list_.begin(); it
				!= action_handler_list_.end(); ++it) {
			applier_->ApplyAdd(*it, edge);
		}
	}

	void FireDeleteVertex(VertexId v) {
		for (auto it = action_handler_list_.rbegin(); it
				!= action_handler_list_.rend(); ++it) {
			applier_->ApplyDelete(*it, v);
		}
	}

	void FireDeleteEdge(EdgeId edge) {
		TRACE("FireDeleteEdge for "<<action_handler_list_.size()<<" handlers");
		for (auto it = action_handler_list_.rbegin(); it
				!= action_handler_list_.rend(); ++it) {
			TRACE("FireDeleteEdge to handler "<<(*it)->name());
			applier_->ApplyDelete(*it, edge);
		}
		TRACE("FireDeleteEdge OK");
	}

	void FireMerge(vector<EdgeId> oldEdges, EdgeId newEdge) {
		TRACE("Fire Merge");
		for (auto it = action_handler_list_.begin(); it
				!= action_handler_list_.end(); ++it) {
			applier_->ApplyMerge(*it, oldEdges, newEdge);
		}
	}

	void FireGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) {
		TRACE("FireGlue for "<<action_handler_list_.size()<<" handlers");
		for (auto it = action_handler_list_.begin(); it
				!= action_handler_list_.end(); ++it) {
			TRACE("FireGlue to handler "<<(*it)->name());
			applier_->ApplyGlue(*it, new_edge, edge1, edge2);
		}
		TRACE("FireGlue OK");
	}

	void FireSplit(EdgeId edge, EdgeId newEdge1, EdgeId newEdge2) {
		TRACE("Fire Split");
		for (auto it = action_handler_list_.begin(); it
				!= action_handler_list_.end(); ++it) {
			applier_->ApplySplit(*it, edge, newEdge1, newEdge2);
		}
	}

public:
	void FireVertexSplit(VertexId newVertex, vector<pair<EdgeId, EdgeId> > newEdges, vector<double> &split_coefficients, VertexId oldVertex) {
		DEBUG("Fire VertexSplit");
		for (auto it = action_handler_list_.begin(); it
				!= action_handler_list_.end(); ++it) {
			applier_->ApplyVertexSplit(*it, newVertex, newEdges, split_coefficients, oldVertex);
		}
	}


	ObservableGraph(HandlerApplier<VertexId, EdgeId> *applier) :
		applier_(applier), element_order_(*this) {
	}

	virtual ~ObservableGraph() {
		TRACE("~ObservableGraph")
		delete applier_;
		TRACE("~ObservableGraph ok")
	}

	void AddActionHandler(Handler* action_handler) const {
		TRACE("Action handler " << action_handler->name() << " added");
		if (find(action_handler_list_.begin(),action_handler_list_.end(), action_handler) != action_handler_list_.end()){
			VERIFY_MSG(false, "Action handler " << action_handler->name() << " has already been added");
		} else {
			action_handler_list_.push_back(action_handler);
		}
	}

	bool RemoveActionHandler(Handler* action_handler) const {
		TRACE("Trying to remove action handler " << action_handler->name());
		for (auto it = action_handler_list_.begin(); it
				!= action_handler_list_.end(); ++it) {
			if (*it == action_handler) {
				action_handler_list_.erase(it);
				TRACE("Action handler " << action_handler->name() << " removed");
				return true;
			}
		}
		TRACE("Action handler " << action_handler->name() << " wasn't found among graph action handlers");
		return false;
	}

	virtual VertexIterator begin() const = 0;

	virtual VertexIterator end() const = 0;

	//todo think of moving to AbstractGraph
	virtual const vector<EdgeId> OutgoingEdges(VertexId vertex) const = 0;

	Comparator ReliableComparatorInstance() const {
		return Comparator(&element_order_);
	}

	template<typename Comparator>
	SmartVertexIterator<ObservableGraph, Comparator> SmartVertexBegin(
			const Comparator& comparator) const {
		return SmartVertexIterator<ObservableGraph, Comparator> (*this,
				comparator);
	}

	SmartVertexIterator<ObservableGraph, Comparator> SmartVertexBegin() const {
		return SmartVertexIterator<ObservableGraph, Comparator> (*this,
				ReliableComparatorInstance());
	}


	template<typename Comparator>
	SmartEdgeIterator<ObservableGraph, Comparator> SmartEdgeBegin(
			const Comparator& comparator) const {
		return SmartEdgeIterator<ObservableGraph, Comparator> (*this,
				comparator);
	}

	SmartEdgeIterator<ObservableGraph, Comparator> SmartEdgeBegin() const {
		return SmartEdgeIterator<ObservableGraph, Comparator> (*this,
				ReliableComparatorInstance());
	}

//	ReliableComparator<VertexId> ReliableComparatorInstance() {
//		return ReliableComparator<VertexId>(element_order_);
//	}


private:
	DECL_LOGGER("ObservableGraph")
};
}
#endif /* OBSERVABLE_GRAPH_HPP_ */
