//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
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

#include "sequential_algorithm.hpp"

#include <memory>

namespace omnigraph {

template <class Graph, class Argument>
class SequentialAlgorihtmFactory {

public:
	typedef std::shared_ptr<SequentialAlgorithm<Argument>> AlgorithmPtr;

	virtual AlgorithmPtr CreateAlgorithm(Graph& graph) = 0;
	virtual ~SequentialAlgorihtmFactory() { }
};

} // namespace omni

#endif /* SEQUENTIAL_ALGORIHTM_FACTORY_HPP_ */
