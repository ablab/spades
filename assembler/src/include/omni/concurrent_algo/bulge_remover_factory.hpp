//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * bulge_remover_factory.hpp
 *
 *  Created on: Sep 15, 2012
 *      Author: alex
 */

#ifndef BULGE_REMOVER_FACTORY_HPP_
#define BULGE_REMOVER_FACTORY_HPP_

#include "omni/sequential_algorihtm_factory.hpp"
#include "omni/concurrent_graph_component.hpp"
#include "omni/sequential_algorihtm_factory.hpp"
#include "omni/bulge_remover.hpp"
#include "kmer_mapper_logger.hpp"

namespace debruijn {


template <class Graph>
class BulgeRemoverFactory
		: public omnigraph::SequentialAlgorihtmFactory<omnigraph::ConcurrentGraphComponent<Graph>, typename Graph::EdgeId> {

	typedef omnigraph::ConcurrentGraphComponent<Graph> Component;
	typedef typename Component::EdgeId EdgeId;
	typedef SequentialAlgorihtmFactory<Component, EdgeId> Base;
	typedef typename Base::AlgorithmPtr AlgorithmPtr;
	typedef typename omnigraph::BulgeRemover<Component>::BulgeCallbackF BulgeCallbackF;
	typedef KmerMapperLogger<Graph> Logger;


public:

	BulgeRemoverFactory(
			size_t max_length,
			double max_coverage,
			double max_relative_coverage,
			double max_delta,
			double max_relative_delta,
			BulgeCallbackF bulge_condition,
			BulgeCallbackF opt_callback = 0,
			boost::function<void(EdgeId)> removal_handler = 0)
				: max_length_(max_length),
				  max_coverage_(max_coverage),
				  max_relative_coverage_(max_relative_coverage),
				  max_delta_(max_delta),
				  max_relative_delta_(max_relative_delta),
				  bulge_condition_(bulge_condition),
				  opt_callback_(opt_callback),
				  removal_handler_(removal_handler) {

		// TODO: KMermapper handling here
	}

	virtual AlgorithmPtr CreateAlgorithm(Component& graph) {
		AlgorithmPtr ptr(
			new BulgeRemover<Component>(
					graph, max_length_, max_coverage_, max_relative_coverage_,
					max_delta_, max_relative_delta_, bulge_condition_,
					opt_callback_, removal_handler_));

		return ptr;
	}

	~BulgeRemoverFactory() {
		while (!loggers_.empty()) {
			delete loggers_.back();
			loggers_.pop_back();
		}
	}

	const vector<Logger*>& loggers() const {
		return loggers_;
	}


private:

	size_t max_length_;
	double max_coverage_;
	double max_relative_coverage_;
	double max_delta_;
	double max_relative_delta_;
	BulgeCallbackF bulge_condition_;
	BulgeCallbackF opt_callback_;
	boost::function<void(EdgeId)> removal_handler_;
	vector<Logger*> loggers_;

};

} // namespace omnigraph

#endif /* BULGE_REMOVER_FACTORY_HPP_ */
