#pragma once

#include "func.hpp"

namespace omnigraph {

template<class Graph>
class ProcessingAlgorithm: private boost::noncopyable {
	Graph& g_;

protected:
	Graph& g() const {
		return g_;
	}

public:
	ProcessingAlgorithm(Graph& g) :
			g_(g) {

	}

	virtual ~ProcessingAlgorithm() {
	}

	virtual bool Process() = 0;
};

template<class Graph, class Comparator = std::less<typename Graph::EdgeId>>
class EdgeProcessingAlgorithm: public ProcessingAlgorithm<Graph> {
	typedef ProcessingAlgorithm<Graph> base;
	typedef typename Graph::EdgeId EdgeId;

	const Comparator comp_;
	const shared_ptr<func::Predicate<EdgeId>> proceed_condition_;

protected:
	virtual bool Process(EdgeId e) = 0;

public:
	EdgeProcessingAlgorithm(Graph& g, const Comparator& c = Comparator(),
			shared_ptr<func::Predicate<EdgeId>> proceed_condition = make_shared<
					func::AlwaysTrue<EdgeId>>()) :
			base(g), comp_(c), proceed_condition_(proceed_condition) {

	}

	bool Process() {
		TRACE("Start processing");
		bool triggered = false;
		for (auto it = this->g().SmartEdgeBegin(comp_); !it.IsEnd(); ++it) {
			if (!proceed_condition_->Check()) {
				TRACE("Stop condition was reached.");
				break;
			}

			TRACE("Processing edge " << this->g().str(*it));
			triggered |= Process(*it);
		}
		TRACE("Finished processing. Triggered = " << triggered);
		return triggered;
	}

private:
	DECL_LOGGER("EdgeProcessingAlgorithm")
	;
};

template<class Graph, class Comparator = std::less<typename Graph::EdgeId>>
class EdgeRemovingAlgorithm: public EdgeProcessingAlgorithm<Graph, Comparator> {
	typedef EdgeProcessingAlgorithm<Graph, Comparator> base;
	typedef typename Graph::EdgeId EdgeId;

	shared_ptr<func::Predicate<EdgeId>> remove_condition_;
	boost::function<void(EdgeId)> removal_handler_;

protected:
	bool Process(EdgeId e) {
		if (remove_condition_->Check(e)) {
			TRACE("Deletion of edge " << this->g().str(e) << " was requested");
			VertexId start = this->g().EdgeStart(e);
			VertexId end = this->g().EdgeEnd(e);

			//todo remove this stupid condition!!!
			if (start == end) {
				return false;
			}

			TRACE("Start " << this->g().str(start));
			TRACE("End " << this->g().str(end));
			if (removal_handler_) {
				TRACE("Calling handler");
				removal_handler_(e);
			}

			TRACE("Deleting edge");
			this->g().DeleteEdge(e);
			TRACE("Compressing locality");
			if (!this->g().RelatedVertices(start, end)) {
				TRACE("Vertices not related");
				TRACE("Compressing end");
				this->g().CompressVertex(end);
				TRACE("End Compressed");
			}
			TRACE("Compressing start");
			this->g().CompressVertex(start);
			TRACE("Start compressed");
			return true;
		} else {
			return false;
		}
	}

public:
	EdgeRemovingAlgorithm(Graph& g,
			shared_ptr<func::Predicate<EdgeId>> remove_condition,
			boost::function<void(EdgeId)> removal_handler = boost::none,
			const Comparator& c = Comparator(),
			shared_ptr<func::Predicate<EdgeId>> proceed_condition = make_shared<
					func::AlwaysTrue<EdgeId>>()) :
			base(g, c, proceed_condition), remove_condition_(remove_condition) {

	}

private:
	DECL_LOGGER("EdgeRemovingAlgorithm")
	;
};

}
