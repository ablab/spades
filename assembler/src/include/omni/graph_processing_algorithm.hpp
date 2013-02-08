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
	virtual bool ProcessEdge(EdgeId e) = 0;

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
			if (!proceed_condition_->Check(*it)) {
				TRACE("Stop condition was reached.");
				break;
			}

			TRACE("Processing edge " << this->g().str(*it));
			triggered |= ProcessEdge(*it);
		}
		TRACE("Finished processing. Triggered = " << triggered);
		return triggered;
	}

private:
	DECL_LOGGER("EdgeProcessingAlgorithm")
	;
};

template<class Graph>
class EdgeRemover {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef boost::function<void(EdgeId)> HandlerF;

	Graph& g_;
	HandlerF removal_handler_;

public:
	EdgeRemover(Graph& g,
			HandlerF removal_handler = 0) :
			g_(g), removal_handler_(removal_handler) {
	}

	bool DeleteEdge(EdgeId e) {
		TRACE("Deletion of edge " << g_.str(e));
		VertexId start = g_.EdgeStart(e);
		VertexId end = g_.EdgeEnd(e);

		TRACE("Start " << g_.str(start));
		TRACE("End " << g_.str(end));
		if (removal_handler_) {
			TRACE("Calling handler");
			removal_handler_(e);
		}
		TRACE("Deleting edge");
		g_.DeleteEdge(e);
		TRACE("Compressing locality");
		if (!g_.RelatedVertices(start, end)) {
			TRACE("Vertices not related");
			TRACE("Compressing end");
			g_.CompressVertex(end);
			TRACE("End Compressed");
		}
		TRACE("Compressing start");
		g_.CompressVertex(start);
		TRACE("Start compressed");
		return true;
	}

private:
	DECL_LOGGER("EdgeRemover")
	;
};

template<class Graph, class Comparator = std::less<typename Graph::EdgeId>>
class EdgeRemovingAlgorithm: public EdgeProcessingAlgorithm<Graph, Comparator> {
	typedef EdgeProcessingAlgorithm<Graph, Comparator> base;
	typedef typename Graph::EdgeId EdgeId;

	shared_ptr<func::Predicate<EdgeId>> remove_condition_;
	EdgeRemover<Graph> edge_remover_;

protected:
	bool ProcessEdge(EdgeId e) {
		if (remove_condition_->Check(e)) {
			return edge_remover_.DeleteEdge(e);
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
			base(g, c, proceed_condition), remove_condition_(remove_condition), edge_remover_(
					g, removal_handler) {

	}

private:
	DECL_LOGGER("EdgeRemovingAlgorithm")
	;
};

}
