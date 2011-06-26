#ifndef OBSERVABLE_GRAPH_HPP_
#define OBSERVABLE_GRAPH_HPP_

#include "omni_utils.hpp"

namespace omnigraph {
template<typename VertexIdT, typename EdgeIdT, typename VertexIterator = typename set<VertexIdT>::iterator>
class ObservableGraph {
public:
	typedef VertexIdT VertexId;
	typedef EdgeIdT EdgeId;
	typedef SmartVertexIterator<ObservableGraph> SmartVertexItarator;
	typedef SmartEdgeIterator<ObservableGraph> SmartEdgeItarator;
private:
	typedef ActionHandler<VertexId, EdgeId> Handler;

	const HandlerApplier<VertexId, EdgeId> *applier_;

	vector<Handler*> action_handler_list_;

protected:
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
		for (auto it = action_handler_list_.begin(); it
				!= action_handler_list_.end(); ++it) {
			applier_->ApplyDelete(*it, v);
		}
	}

	void FireDeleteEdge(EdgeId edge) {
		TRACE("FireDeleteEdge for "<<action_handler_list_.size()<<" handlers");
		for (auto it = action_handler_list_.begin(); it
				!= action_handler_list_.end(); ++it) {
			TRACE("FireDeleteEdge to handler "<<*it);
			applier_->ApplyDelete(*it, edge);
		}
		TRACE("FireDeleteEdge OK");
	}

	void FireMerge(vector<EdgeId> oldEdges, EdgeId newEdge) {
		for (auto it = action_handler_list_.begin(); it
				!= action_handler_list_.end(); ++it) {
			applier_->ApplyMerge(*it, oldEdges, newEdge);
		}
	}

	void FireGlue(EdgeId edge1, EdgeId edge2) {
		for (auto it = action_handler_list_.begin(); it
				!= action_handler_list_.end(); ++it) {
			applier_->ApplyGlue(*it, edge1, edge2);
		}
	}

	void FireSplit(EdgeId edge, EdgeId newEdge1, EdgeId newEdge2) {
		for (auto it = action_handler_list_.begin(); it
				!= action_handler_list_.end(); ++it) {
			applier_->ApplySplit(*it, edge, newEdge1, newEdge2);
		}
	}


public:

	ObservableGraph(HandlerApplier<VertexId, EdgeId> *applier) :
		applier_(applier) {
	}

	virtual ~ObservableGraph() {
		TRACE("~ObservableGraph")
		delete applier_;
		TRACE("~ObservableGraph ok")
	}

	void AddActionHandler(Handler* action_handler) {
		TRACE("Action handler added");
		if (find(action_handler_list_.begin(),action_handler_list_.end(), action_handler) != action_handler_list_.end()){
			TRACE("Action handler already presented");
//			assert(0);
		}
		else
		action_handler_list_.push_back(action_handler);
	}

	bool RemoveActionHandler(Handler* action_handler) {
		TRACE("Trying to remove action handler");
		for (auto it = action_handler_list_.begin(); it
				!= action_handler_list_.end(); ++it) {
			if (*it == action_handler) {
				action_handler_list_.erase(it);
				TRACE("Action handler removed. Remain "<<action_handler_list_.size());
				return true;
			}
		}
		return false;
	}

	virtual VertexIterator begin() const = 0;

	virtual VertexIterator end() const = 0;

	virtual vector<EdgeId> OutgoingEdges(VertexId vertex) const = 0;

	template<typename Comparator = std::less<VertexId> >
	SmartVertexIterator<ObservableGraph, Comparator> SmartVertexBegin(
			const Comparator& comparator = Comparator()) {
		return SmartVertexIterator<ObservableGraph, Comparator> (*this,
				true, comparator);
	}

	template<typename Comparator = std::less<VertexId> >
	SmartVertexIterator<ObservableGraph, Comparator> SmartVertexEnd(
			const Comparator& comparator = Comparator()) {
		return SmartVertexIterator<ObservableGraph, Comparator> (*this,
				false, comparator);
	}

	template<typename Comparator = std::less<EdgeId> >
	SmartEdgeIterator<ObservableGraph, Comparator> SmartEdgeBegin(
			const Comparator& comparator = Comparator()) {
		return SmartEdgeIterator<ObservableGraph, Comparator> (*this,
				true, comparator);
	}

	template<typename Comparator = std::less<EdgeId> >
	SmartEdgeIterator<ObservableGraph, Comparator> SmartEdgeEnd(
			const Comparator& comparator = Comparator()) {
		return SmartEdgeIterator<ObservableGraph, Comparator> (*this,
				false, comparator);
	}


private:
	DECL_LOGGER("ObservableGraph")
};

}
#endif /* OBSERVABLE_GRAPH_HPP_ */
