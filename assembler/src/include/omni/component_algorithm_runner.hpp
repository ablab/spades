//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
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

#include <boost/shared_ptr.hpp>

#include "concurrent_graph_component.hpp"
#include "sequential_algorithm.hpp"

namespace omnigraph {

template <class Graph, class Argument>
class ComponentAlgorithmRunner {

public:
	typedef ConcurrentGraphComponent<Graph> Component;
	typedef boost::shared_ptr<SequentialAlgorithm<Argument>> AlgorithmPtr;


	ComponentAlgorithmRunner(Component& component, AlgorithmPtr algorithm)
			: component_(component),
			  algorithm_(algorithm),
			  not_processed_arguments_(component, std::less<Argument>(), false) {
	}

	template <class JavaStyleIterator>
	void Run(JavaStyleIterator& it) {
		algorithm_->Preprocessing();

		VERIFY(component_.IsValid());

		for (; !it.IsEnd(); ++it) {
			algorithm_->ProcessNext(*it);

			if (!component_.IsValid()) {
				not_processed_arguments_.insert(*it);
				component_.ValidateComponent();
			}
		}

		algorithm_->Postprocessing();

		VERIFY(component_.IsValid());
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
