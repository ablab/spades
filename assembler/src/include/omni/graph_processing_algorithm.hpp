#pragma once

#include "func.hpp"

namespace omnigraph {

template <class Graph>
class ProcessingAlgorithm : private boost::noncopyable {
	Graph& g_;

protected:
	Graph& g() const {
		return g_;
	}

public:
	ProcessingAlgorithm(Graph& g) : g_(g) {

	}

	virtual ~ProcessingAlgorithm() { }

	virtual bool Process() = 0;
};

template <class Graph, class Comparator = std::less<typename Graph::EdgeId>>
class EdgeProcessingAlgorithm : public ProcessingAlgorithm<Graph> {
	typedef ProcessingAlgorithm<Graph> base;
	typedef typename Graph::EdgeId EdgeId;

	const Comparator comp_;
	const shared_ptr<func::Predicate<EdgeId>> proceed_condition_;

protected:
	virtual bool Process(EdgeId e) = 0;

public:
	EdgeProcessingAlgorithm(Graph& g, const Comparator& c, shared_ptr<func::Predicate<EdgeId>> proceed_condition
			= make_shared<func::AlwaysTrue<EdgeId>>()) : base(g), comp_(c), proceed_condition_(proceed_condition) {

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
			triggered = Process(*it);
		}
		TRACE("Finished processing. Triggered = " << triggered);
		return triggered;
	}

private:
	DECL_LOGGER("EdgeProcessingAlgorithm");
};

template <class Graph, class Comparator = std::less<typename Graph::EdgeId>>
class EdgeRemovingAlgorithm : public EdgeProcessingAlgorithm {
	typedef EdgeProcessingAlgorithm<Graph> base;
	typedef typename Graph::EdgeId EdgeId;

protected:
	bool Process(EdgeId e) {

	}

public:
	EdgeRemovingAlgorithm(Graph& g, const Comparator& c, shared_ptr<func::Predicate<EdgeId>> remove_condition
			, shared_ptr<func::Predicate<EdgeId>> proceed_condition)
	: base(g, c, proceed_condition) {

	}

private:
	DECL_LOGGER("EdgeRemovingAlgorithm");
};

}
