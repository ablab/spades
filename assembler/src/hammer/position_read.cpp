//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "standard.hpp"
#include "position_kmer.hpp"
#include "position_read.hpp"

char PositionRead::at(size_t pos) const {
	return Globals::blob[start_ + pos];
}

char PositionRead::operator[](size_t pos) const {
	return Globals::blob[start_ + pos];
}

pair<size_t, hint_t> PositionRead::nextKMerNo(size_t begin) const {
	for (size_t pos = begin + 1; pos < size_ - K + 1; ++pos) {
		//cout << "    looking for " << (start_ + pos) << ": " << PositionKMer(start_ + pos).str();
		vector<hint_t>::const_iterator it_vec = lower_bound(Globals::kmernos->begin(), Globals::kmernos->end(), start_ + pos, PositionKMer::compareKMersDirect);
		//cout << "  lowerbound=" << (it_vec - Globals::kmernos->begin()) << "\tvalue=" << *it_vec;
		if (it_vec != Globals::kmernos->end() && PositionKMer::equalKMersDirect(*it_vec, start_ + pos)) {
			//cout << "  ok" << endl;
			return make_pair(pos, (it_vec - Globals::kmernos->begin()));
		}
	}
	//cout << endl;

	return make_pair(-1, BLOBKMER_UNDEFINED);
}
