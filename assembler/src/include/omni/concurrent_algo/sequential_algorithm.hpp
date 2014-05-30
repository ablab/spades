//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************


/*
 * sequential_algorithm.hpp
 *
 *  Created on: Sep 7, 2012
 *      Author: Alexander Opeykin (alexander.opeykin@gmail.com)
 */

/*
 * Interface for algorithms that can be processed part by part.
 */

#ifndef SEQUENTIAL_ALGORITHM_HPP_
#define SEQUENTIAL_ALGORITHM_HPP_

namespace omnigraph {

template <class T>
class SequentialAlgorithm {

public:
	virtual ~SequentialAlgorithm() { }

	virtual void Preprocessing() { }
	virtual void Postprocessing() { }
	virtual bool ProcessNext(const T& arg) = 0;
};

} //namespace omnigraph
#endif /* SEQUENTIAL_ALGORITHM_HPP_ */
