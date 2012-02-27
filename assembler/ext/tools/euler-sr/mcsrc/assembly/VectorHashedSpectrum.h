/***************************************************************************
 * Title:          HashedSpectrum.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  02/27/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef HASHED_SPECTRUM_H_
#define HASHED_SPECTRUM_H_

#include "Spectrum.h"
#include "ReadPos.h"
#include "SimpleSequence.h"
#include "IntegralTuple.h"
#include "hash/VectorHashTable.h" // this defines the storage



template <typename T_TupleType, ssize_t T_HashLength, ssize_t BUFFER_SIZE>
class VectorTupleHash  {
 public:
	// use for computing a hash
	IntegralTupleHashValueFunctor<T_HashLength>  hashFunction;
	UpdateIntegralTupleCountFunctor<T_TupleType> updateFunction;

	// store here
	typedef VectorHashTable<T_TupleType,
		IntegralTupleHashValueFunctor<T_HashLength>,
		UpdateIntegralTupleCountFunctor<T_TupleType>, BUFFER_SIZE > HashTable;

	HashTable hashTable;

	void Flush() {
		hashTable.Flush(hashFunction, updateFunction);
	}

	void FlushDirected() {
		hashTable.FlushDirected(hashFunction, updateFunction);
	}

	VectorTupleHash() {
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
		ssize_t tupleSet = 0;
		for (p = trimFront; p < seq.length - T_TupleType::tupleSize + 1 - trimEnd; p++) {
			if (unmasked_nuc_index[seq.seq[p + T_TupleType::tupleSize -1]] >= 4) {
				nextInvalid = p + T_TupleType::tupleSize -1;
				tupleSet = 0;
			}
			if (p > nextInvalid) {
				if (tupleSet == 0) {
					tuple.StringToTuple(&(seq.seq[p]));
					tupleSet = 1;
				}
				else {
					T_TupleType nextTuple;
					ForwardNuc(tuple, 
										 unmasked_nuc_index[seq.seq[p+IntegralTuple::tupleSize - 1]], 
										 nextTuple);
					//					tuple.tuple = nextTuple.tuple;
					tuple.CopyTuple(nextTuple);
				}
				if (hashTable.BufferedStore(tuple, hashFunction, updateFunction))
					newTupleStored = 1;
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
		//UNUSED+// T_TupleType* *revStrand;
		T_TupleType* forStrand ;
		forStrand = 0;
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
				hashTable.BufferedStore(tuple, hashFunction, updateFunction);

				tupleCopy = tuple;
				//				std::cout << "original tuple: " << tuple.tuple << std::endl;
				query = tuple;
				tuple.MakeRC(queryRC);
				queryRCCopy = queryRC;
				forStrand = hashTable.Find(query, hashFunction);
				//				revStrand = hashTable.Find(queryRC, hashFunction);
				//				std::cout << query.tuple << " " << forStrand << " " << revStrand << std::endl;
				if (forStrand) {
					// The key already exists in the forward direction, update that one
					updateFunction(*forStrand);
					//					if (hashTable.BufferedStore(tuple, hashFunction, updateFunction))
					//					if (hashTable.Store(tuple, hashFunction, updateFunction))
					newTupleStored = 1;
				}
			}
		}
		return newTupleStored;
	}
};

#endif
