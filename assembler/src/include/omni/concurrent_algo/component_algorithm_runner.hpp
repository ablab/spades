//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************


/*
 * component_edge_algorithm.hpp
 *
 *  Created on: Sep 7, 2012
 *      Author: Alexander Opeykin (alexander.opeykin@gmail.com)
 */


#ifndef COMPONENT_ALGORITHM_RUNNER_HPP_
#define COMPONENT_ALGORITHM_RUNNER_HPP_

#include <memory>

#include "concurrent_graph_component.hpp"
#include "sequential_algorithm.hpp"

namespace omnigraph {

//Deprecated! Use SmartSetIterator instead!
template<class Graph, typename ElementId, typename Comparator = std::less<
        ElementId> >
class SmartSet: public GraphActionHandler<Graph> {
public:
    typedef typename set<ElementId, Comparator>::iterator iterator;
    typedef typename set<ElementId, Comparator>::const_iterator const_iterator;
private:
    set<ElementId, Comparator> inner_set_;
    const bool add_new_;

public:
    SmartSet(const Graph &graph, Comparator comparator = Comparator(),
            bool add_new = true) :
            GraphActionHandler<Graph>(graph, "SmartSet"), inner_set_(
                    comparator), add_new_(add_new) {
    }

    template<class Iter>
    SmartSet(Iter begin, Iter end, const Graph &graph, Comparator comparator =
            Comparator(), bool add_new = true) :
            GraphActionHandler<Graph>(graph, "SmartSet"), inner_set_(begin, end,
                    comparator), add_new_(add_new) {
    }

    virtual ~SmartSet() {
    }

    virtual void HandleAdd(ElementId v) {
        if (add_new_)
            inner_set_.insert(v);
    }

    virtual void HandleDelete(ElementId v) {
        inner_set_.erase(v);
    }

    iterator begin() {
        return inner_set_.begin();
    }

    iterator end() {
        return inner_set_.end();
    }

    const_iterator begin() const {
        return inner_set_.begin();
    }

    const_iterator end() const {
        return inner_set_.end();
    }

    pair<iterator, bool> insert(const ElementId& elem) {
        return inner_set_.insert(elem);
    }

    const set<ElementId, Comparator> &inner_set() {
        return inner_set_;
    }
};

template <class Graph, class Argument>
class ComponentAlgorithmRunner {

public:
	typedef ConcurrentGraphComponent<Graph> Component;
	typedef std::shared_ptr<SequentialAlgorithm<Argument>> AlgorithmPtr;


	ComponentAlgorithmRunner(Component& component, AlgorithmPtr algorithm)
			: component_(component),
			  algorithm_(algorithm),
			  not_processed_arguments_(component, std::less<Argument>(), false) {
	}

	template <class JavaStyleIterator>
	void Run(JavaStyleIterator& it) {
		algorithm_->Preprocessing();

		for (; !it.IsEnd(); ++it) {
			if (!algorithm_->ProcessNext(*it)) {
				not_processed_arguments_.insert(*it);
			}
		}

		algorithm_->Postprocessing();
	}

	void GetNotProcessedArguments(vector<Argument>& arguments) {
		arguments.insert(arguments.end(), not_processed_arguments_.begin(), not_processed_arguments_.end());
	}

private:
	Component& component_;
	AlgorithmPtr algorithm_;
	SmartSet<Component, Argument> not_processed_arguments_;
};

} // namespace omnigraph



#endif /* COMPONENT_ALGORITHM_HPP_ */
