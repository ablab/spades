/***************************************************************************
 * Title:          VectorHashTable.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  02/27/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef VECTOR_HASH_TABLE
#define VECTOR_HASH_TABLE

#include <vector>



template<typename T, typename HashFunct, typename UpdateFunct, ssize_t BUFFER_SIZE>
class VectorHashTable {
public:
  ssize_t size;
	ssize_t count;
  vector<T> **table;
	T buffer[BUFFER_SIZE];
	ssize_t bufCur;

	VectorHashTable() { count = 0;};
  VectorHashTable(HashFunct &hashFunct) {
		count = 0;
		Init(hashFunct);
  }
	void Init(HashFunct &hashFunct) {
    size = hashFunct.MaxHashValue();
    table = new vector<T>*[size];
		fill(table, table + size, (vector<T>*) NULL);
		bufCur = 0;
	}

  ssize_t CountSize() {
    ssize_t i;
    ssize_t total = 0;
    for (i = 0; i < size; i++ ) {
			if (table[i] != NULL) {
				total += table[i]->size();
			}
		}
    return total;
  }
	class SortOnHash {
	public:
		HashFunct hash;
		ssize_t operator()(const T a, const T b) const {
			return this->hash(a) < this->hash(b);
		}
	};

#if 0 // Never called, so disable
	class SortOnNHash {
	public:
		ssize_t hashLength;
		LongTuple mask;
		HashFunct & hashFunct;
		ssize_t operator()(const LongTuple a, const LongTuple b) const {
			// TODO: Do we need to use comparison for the new LongTuple?
			// It appears this is never called, so it doesn't matter.
			return (hashFunct(a) < hashFunct(b));
		}
	};
#endif
			
	ssize_t FlushDirected(HashFunct &hashFunct, UpdateFunct &update) {
		//		cout << "Flushing " << bufCur << endl;
		ssize_t i = 0;
		SortOnHash sortOnHash;
		sortOnHash.hash = hashFunct;
		std::sort(&buffer[0], &buffer[bufCur], sortOnHash);
		//		std::sort(&buffer[0], &buffer[bufCur]);
		// Compress the buffer
		if (bufCur == 0)
			return 1;
		/*
		 * Compute some performance diagnostics.
		ssize_t numCacheHits = 0;
		ssize_t curIndex = 0;
		for (i = 0; i < bufCur; i++) {
			curIndex = i;
			while (i < bufCur and hashFunct(buffer[curIndex]) == hashFunct(buffer[i].tuple)) {
				i++;
			}
			if (i > curIndex) {
				numCacheHits++;
			}
		}		
		cout << "num cache hits: " << numCacheHits << std::endl;

		ssize_t numCompressed = 0;
		T cur;
		ssize_t curIndex = 0;
		for (i = 0; i < bufCur; i++) {
			cur = buffer[i];
			curIndex = i;
			while (i < bufCur and cur.tuple == buffer[i].tuple) {
				update(cur);
				i++;
			}
			buffer[curIndex] = cur;
			curIndex++;
			if (curIndex > i) numCompressed++;
		}
		bufCur = curIndex;
		cout << "Compressed " << numCompressed << " / " << bufCur << endl;
		*/
		for (i = 0; i < bufCur; i++) {
			T* queryLookupAddr;
			if ((queryLookupAddr = Find(buffer[i])) != NULL) {
				update(*queryLookupAddr);
				//				Store(buffer[i], hashFunct, update);
			}
			else {
				T queryRC, queryRCCopy;
				buffer[i].MakeRC(queryRC);
				queryRCCopy = queryRC;
				if ((queryLookupAddr = Find(queryRC)) != NULL)
					update(*queryLookupAddr);
				else
					StoreNew(buffer[i], hashFunct);
			}
		}
		bufCur = 0;
		return 1;
	}

	ssize_t Flush(HashFunct &hashFunct, UpdateFunct &update) {
		//		cout << "Flushing " << bufCur << endl;
		ssize_t i = 0;
		std::sort(&buffer[0], &buffer[bufCur]);
		// Compress the buffer
		if (bufCur == 0)
			return 1;
		T cur;
		cur = buffer[0];
		ssize_t numCompressed = 0;
		for (i = 1; i < bufCur; i++) {
			if (IntegralTuple::cmp(cur.tuple,buffer[i].tuple) == 0) {
				numCompressed++;
				while (i < bufCur and IntegralTuple::cmp(cur.tuple,buffer[i].tuple) == 0)
					i++;
			}
#if 0 // OLD CODE
			if (cur.tuple == buffer[i].tuple) {
				numCompressed++;
				while (i < bufCur and cur.tuple == buffer[i].tuple)
					i++;
			}
#endif
			cur = buffer[i];
		}
		cout << "Compressed " << numCompressed << " / " << bufCur << endl;
		for (i = 0; i < bufCur; i++) {
			Store(buffer[i], hashFunct, update);
		}
		bufCur = 0;
		return 1;
	}

	ssize_t BufferedStore(T &value, HashFunct &hashFunct, UpdateFunct &update) {
		buffer[bufCur] = value;
		bufCur++;
		if (bufCur >= BUFFER_SIZE) {
			FlushDirected(hashFunct, update);
		}
		return 1; // TODO: return value is being checked for success; but when can this fail?
	}

	ssize_t StoreNew(T &value, HashFunct &hashFunct) {
    ssize_t index = hashFunct(value);
    if (index < 0) 
     return 0;
		
		// If this is new, simply return here.
		if (table[index] == NULL) {
			table[index] = new vector<T>;
		}
		
		table[index]->push_back(value);
		++count;
		return 1;
	}
		
		
  ssize_t Store(T &value, HashFunct &hashFunct, UpdateFunct &update) {
    ssize_t index = hashFunct(value);
    if (index < 0) 
     return 0;
		
		// If this is new, simply return here.
		if (table[index] == NULL) {
			table[index] = new vector<T>;
			table[index]->push_back(value);
			++count;
			return 1;
		}

		// The table must have a vector here
		assert(table[index] != NULL);
		typename vector<T>::iterator it;
		it = std::find(table[index]->begin(), table[index]->end(), value);
		if (it != table[index]->end() and *it == value) {
			update(*it);
			return 0;
		}
		else {
			table[index]->push_back(value);
			//			std::sort(table[index]->begin(), table[index]->end());
			count++;
			return 1;
		}
  }
	
	T* Find(T &value, HashFunct hashFunct=  HashFunct()) {
		
		ssize_t index = hashFunct(value);
		assert(index < size);
		//		cout << "size: " << size << " index: " << index << endl;
		if (table[index] == NULL) {
			return NULL;
		}
		else {
			typename vector<T>::iterator it;
			it = std::find(table[index]->begin(), table[index]->end(), value);
			if (it != table[index]->end() and *it == value) {
				value = *it;
				return (T*) &*it;//table[index]->begin() + (it - table[index]->begin()); // strange math needed to go back from iterator to addresss
			}
			else {
				return NULL;
			}
		}
	}
		
	void Summarize() {
		ssize_t nChains = 0;
		ssize_t chainLength, totalChainLength = 0;
		ssize_t i;
		ssize_t maxLen = 200;
		ssize_t bins[maxLen];
		ssize_t b;

		for (b = 0; b < maxLen; b++ ) bins[b] = 0;

		for (i =0 ; i < size; i++) {
			if (table[i] != NULL) {
				nChains++;
				chainLength = table[i]->size();
				if (chainLength  > maxLen-1) 
					bins[maxLen-1]++;
				else 
					bins[chainLength]++;
				totalChainLength += table[i]->size();
			}
		}
		for (b = 0; b < maxLen; b++ )
			std::cout << bins[b] << " ";
		cout << endl;
		std::cout << nChains << " vectors." << std::endl;
		std::cout << totalChainLength / nChains << " average vector length. " << std::endl;
	}

	void Free() {
		ssize_t i;
		for (i = 0; i < size; i++) {
			if (table[i] != NULL)
				delete table[i];
		}
	}
};
    

#endif
