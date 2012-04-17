//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * read_reconstruction.hpp
 *
 *  Created on: Nov 28, 2011
 *      Author: snikolenko
 */

#ifndef READ_RECONSTRUCTION_HPP_
#define READ_RECONSTRUCTION_HPP_

#include "globals.hpp"

class ReadReconstructor {
public:
	ReadReconstructor() : read_corrected(Globals::pr->size(), false), changed_reads(0), changed_bases(0) {
	}

	void reconstruct(const std::string & fname);
	void reconstruct(const std::string & fn_left, const std::string & fn_right);

	hint_t getChangedReads() { return changed_reads; }
	hint_t getChangedBases() { return changed_bases; }

private:
	std::vector<Read> reads;
	std::vector<bool> read_corrected;
	hint_t changed_reads;
	hint_t changed_bases;

	void reconstructAll();
	bool makeOneReconstructionIteration();
	bool correctOneRead(PositionRead &);
};

#endif /* READ_RECONSTRUCTION_HPP_ */
