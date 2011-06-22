#ifndef OBSERVABLE_GRAPH_HPP_
#define OBSERVABLE_GRAPH_HPP_

#include "omni_utils.hpp"

namespace omnigraph {
template<typename VertexId, typename EdgeId>
class ObservableGraph {
private:
	typedef ActionHandler<VertexId, EdgeId> Handler;

	const HandlerApplier<VertexId, EdgeId> *applier_;

	vector<Handler*> action_handler_list_;

protected:
	void FireAddVertex(VertexId v) {
		for (auto it = action_handler_list_.begin(); it
				!= action_handler_list_.end(); ++it) {
			applier_.ApplyAdd(*it, v);
		}
	}

	void FireAddEdge(EdgeId edge) {
		for (auto it = action_handler_list_.begin(); it
				!= action_handler_list_.end(); ++it) {
			applier_.ApplyAdd(*it, edge);
		}
	}

	void FireDeleteVertex(VertexId v) {
		for (auto it = action_handler_list_.begin(); it
				!= action_handler_list_.end(); ++it) {
			applier_->ApplyDelete(*it, v);
		}
	}

	void FireDeleteEdge(EdgeId edge) {
		for (auto it = action_handler_list_.begin(); it
				!= action_handler_list_.end(); ++it) {
			applier_->ApplyDelete(*it, edge);
		}
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

	ObservableGraph(HandlerApplier<VertexId, EdgeId> *applier) :
		applier_(applier) {
	}

	virtual ~ObservableGraph() {
		delete applier_;
	}

	void AddActionHandler(Handler* action_handler) {
		TRACE("Action handler added");
		action_handler_list_.push_back(action_handler);
	}

	bool RemoveActionHandler(Handler* action_handler) {
		TRACE("Trying to remove action handler");
		for (auto it = action_handler_list_.begin(); it
				!= action_handler_list_.end(); ++it) {
			if (*it == action_handler) {
				action_handler_list_.erase(it);
				TRACE("Action handler removed");
				return true;
			}
		}
		return false;
	}

private:
	DECL_LOGGER("ObservableGraph")
};

}
#endif /* OBSERVABLE_GRAPH_HPP_ */
