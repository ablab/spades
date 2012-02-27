/***************************************************************************
 * Title:          StringListSpectrum.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef LIST_SPECTRUM_H_
#define LIST_SPECTRUM_H_

#include "Spectrum.h"

template<typename T>
class ListSpectrum : public Spectrum<T> {
public:
  std::vector<T> tupleList;

	ListSpectrum() {
		this->tupleSize = -1;
	}
	void Resize(ssize_t size) {
		tupleList.resize(size);
	}
	ssize_t size() {
		return tupleList.size();
	}
	ssize_t ReadSpectrum(std::string &fileName, ssize_t minMult = 0);
	void WriteSpectrum(std::string &fileName, ssize_t minMult = 0);
	ssize_t IncrementKmer(T &tuple);
	ssize_t FindKmer(Tuple &tuple);
	ssize_t FindKmer(char *valueStr);
};



#endif

