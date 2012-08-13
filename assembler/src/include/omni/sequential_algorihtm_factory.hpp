//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************


/*
 * sequential_algorihtm_factory.hpp
 *
 *  Created on: Sep 7, 2012
 *      Author: Alexander Opeykin (alexander.opeykin@gmail.com)
 */


#ifndef SEQUENTIAL_ALGORIHTM_FACTORY_HPP_
#define SEQUENTIAL_ALGORIHTM_FACTORY_HPP_

#include <boost/shared_ptr.hpp>

#include "sequential_algorithm.hpp"

namespace omnigraph {

template <class Graph, class Argument>
class SequentialAlgorihtmFactory {

public:
	typedef boost::shared_ptr<SequentialAlgorithm<Argument>> AlgorithmPtr;

	virtual AlgorithmPtr CreateAlgorithm(Graph& graph) = 0;
	virtual ~SequentialAlgorihtmFactory() { }
};

} // namespace omni

#endif /* SEQUENTIAL_ALGORIHTM_FACTORY_HPP_ */
