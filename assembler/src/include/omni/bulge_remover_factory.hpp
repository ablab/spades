/*
 * bulge_remover_factory.hpp
 *
 *  Created on: Sep 15, 2012
 *      Author: alex
 */

#ifndef BULGE_REMOVER_FACTORY_HPP_
#define BULGE_REMOVER_FACTORY_HPP_

#include "sequential_algorihtm_factory.hpp"
#include "bulge_remover.hpp"

namespace omnigraph {


template <class Graph>
class BulgeRemoverFactory : public SequentialAlgorihtmFactory<Graph, typename Graph::EdgeId> {

	typedef typename Graph::EdgeId EdgeId;
	typedef SequentialAlgorihtmFactory<Graph, EdgeId> Base;
	typedef typename Base::AlgorithmPtr AlgorithmPtr;
	typedef typename omnigraph::BulgeRemover<Graph>::BulgeCallbackF BulgeCallbackF;


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
	}

	virtual AlgorithmPtr CreateAlgorithm(Graph& graph) {
		AlgorithmPtr ptr(
			new BulgeRemover<Graph>(
					graph, max_length_, max_coverage_, max_relative_coverage_,
					max_delta_, max_relative_delta_, bulge_condition_,
					opt_callback_, removal_handler_));

		return ptr;
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

};

} // namespace omnigraph

#endif /* BULGE_REMOVER_FACTORY_HPP_ */
