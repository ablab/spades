/***************************************************************************
 * Title:          HashedSpectrum.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  02/18/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef HASHED_SPECTRUM_H_
#define HASHED_SPECTRUM_H_

#include "Spectrum.h"
#include "ReadPos.h"
#include "SimpleSequence.h"
#include "PagedList.h"
#include "hash/PagedHashTable.h"
#include "IntegralTuple.h"

template <typename T_TupleType, ssize_t T_HashLength, ssize_t T_PageSize>
class TupleHash  {
 public:
	// use for computing a hash
	IntegralTupleHashValueFunctor<T_HashLength>  hashFunction;
	UpdateIntegralTupleCountFunctor<T_TupleType> updateFunction;

	// store here
	typedef PagedHashTable<T_TupleType, T_PageSize, 
		IntegralTupleHashValueFunctor<T_HashLength>,
		UpdateIntegralTupleCountFunctor<T_TupleType> > HashTable;

	HashTable hashTable;

	
	TupleHash() {
		hashTable.Init(hashFunction);
	}

	ssize_t HashSequence(SimpleSequence &seq,
									 ssize_t trimFront = 0,
									 ssize_t trimEnd   = 0) {
		ssize_t p;
		ssize_t newTupleStored = 0;
		ssize_t nextInvalid = -1;
		// Find the first position that is not valid
		for (p = trimFront; p < seq.length - T_TupleType::tupleSize; p++) {
			if (unmasked_nuc_index[seq.seq[p]] >= 4)
				nextInvalid = p;
		}
		T_TupleType tuple;
		for (p = trimFront; p < seq.length - T_TupleType::tupleSize + 1 - trimEnd; p++) {
			if (unmasked_nuc_index[seq.seq[p + T_TupleType::tupleSize -1]] >= 4)
				nextInvalid = p + T_TupleType::tupleSize -1;
			if (p > nextInvalid) {
				tuple.StringToTuple(&(seq.seq[p]));
				if (hashTable.Find(tuple, hashFunction)) {
					if (hashTable.Store(tuple, hashFunction, updateFunction))
						newTupleStored = 1;
				}
			}
		}
		return newTupleStored;
	}

	ssize_t HashSequenceUndirected(SimpleSequence &seq,
														 ssize_t trimFront = 0,
														 ssize_t trimEnd   = 0) {
		ssize_t p;
		ssize_t newTupleStored = 0;
		ssize_t nextInvalid = -1;
		ssize_t forStrand, revStrand;
		forStrand = 0; revStrand = 0;
		// Find the first position that is not valid
		for (p = trimFront; p < seq.length - T_TupleType::tupleSize; p++) {
			if (unmasked_nuc_index[seq.seq[p]] >= 4)
				nextInvalid = p;
		}
		T_TupleType tuple, query, queryRC, tupleCopy, queryRCCopy;
		for (p = trimFront; p < seq.length - T_TupleType::tupleSize + 1 - trimEnd; p++) {
			
			if (unmasked_nuc_index[seq.seq[p + T_TupleType::tupleSize -1]] >= 4)
				nextInvalid = p + T_TupleType::tupleSize -1;
			if (p > nextInvalid) {
				tuple.StringToTuple(&(seq.seq[p]));
				tupleCopy = tuple;
				//				std::cout << "original tuple: " << tuple.tuple << std::endl;
				query = tuple;
				tuple.MakeRC(queryRC);
				queryRCCopy = queryRC;
				forStrand = hashTable.Find(query, hashFunction);
				revStrand = hashTable.Find(queryRC, hashFunction);
				//				std::cout << query.tuple << " " << forStrand << " " << revStrand << std::endl;
				if (forStrand) {
					if (hashTable.Store(tuple, hashFunction, updateFunction))
						newTupleStored = 1;
				}
				else if (revStrand) {
					if (hashTable.Store(queryRCCopy, hashFunction, updateFunction))
						newTupleStored = 1;
				}
				else {
					hashTable.Store(tuple, hashFunction, updateFunction);
				}
			}
		}
		return newTupleStored;
	}
};

#endif
