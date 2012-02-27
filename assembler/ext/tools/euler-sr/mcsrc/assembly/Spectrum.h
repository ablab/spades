/***************************************************************************
 * Title:          Spectrum.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef SPECTRUM_H_
#define SPECTRUM_H_

#include <string>
#include "DNASequence.h"
#include "DeBruijnGraph.h"
#include "ReadPos.h"

// Define the interface for a spectrum
template <typename T>
class Spectrum {
public:
	typedef T TupleType;
  int tupleSize;
	virtual ssize_t size()= 0;
	virtual void Read(std::string &fileName, ssize_t minMult = 0) = 0;
	virtual void Write(std::string &fileName, ssize_t minMult = 0)= 0;
	virtual ssize_t IncrementMult(T &tuple) = 0;
	virtual ssize_t FindTuple(T &tuple) = 0;
};





//  Look for value in 'sortedKmers'.  If value is present, return
//  its index in sortedKmers.  Otherwise, return -1;





#endif
