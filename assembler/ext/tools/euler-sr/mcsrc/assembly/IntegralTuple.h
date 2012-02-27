/***************************************************************************
 * Title:          IntegralTuple.h 
 * Author:         Mark Chaisson, Glenn Tesler
 * Created:        2007
 * Last modified:  03/02/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

// For using in libraries, or when included from other .h files:
//     #include "IntegralTuple.h"
// For using in executables:
//     #include "IntegralTupleStatic.h"

#ifndef INTEGRAL_TUPLE_H_
#define INTEGRAL_TUPLE_H_

#include <iostream>
#include <string>
#include <sstream>

#include "SeqUtils.h"
#include "utils.h"
#include "compatibility.h"

#define MAX_COUNT 1000000

#define LT_BitsPerWord   _NBITS_T_(LongTupleWord)
#define LT_NucsPerWord   ( LT_BitsPerWord / 2 )

// LongTuple is an array of LongTupleWord
// nucleotide 0: word 0, bits 0..1
// nucleotide 1: word 1, bits 2..3
// ...
// LT_WordNum(i) returns the word # with nucleotide i
// LT_BitNum(i) & LT_BitNum(i)+1 are the bit numbers in that word
#define LT_WordNum(i)    ( (i) / LT_NucsPerWord )
#define LT_BitNum(i)     ( 2*((i) % LT_NucsPerWord) )

// LT_Mask(i,nuc)
//   where i has already been reduced to within one word 0..(LT_NucsPerWord-1)
//   This gives nuc shifted to the ith nucleotide in the word.
// LT_Mask_Safe(i,nuc)
//   Returns 0 when i >= LT_NucsPerWord.
// LT_Mask_Mod(i,nuc)
//   Mask for nuc, shifted by (i modulo LT_NucsPerWord)
#define LT_1             LongTupleWord(1)
#define LT_3             LongTupleWord(3)
//#define LT_Mask1(i)      ( LT_1 << (2*(i)) )
//#define LT_Mask3(i)      ( LT_3 << (2*(i)) )
#define LT_Mask(i,nuc)   ( LongTupleWord(nuc) << (2*(i)) )
#define LT_Mask_Mod(i,nuc)   ( LongTupleWord(nuc) << (2*((i) % LT_NucsPerWord)) )
#define LT_Mask_Safe(i,nuc)   ( ((i) < LT_NucsPerWord) \
																? (LongTupleWord(nuc) << (2*(i))) \
																: LongTupleWord(0) )

using namespace std;

class IntegralTuple  {
 public:
	// tuple location
	static int tupleSize;            // number of nucleotides
	static int numWords;             // number of LongTupleWord words needed
	static int numWords1;            // numWords - 1 (index of high order word)
	//	static LongTuple MASK_OFF_FIRST; // no longer used
	// REDO????
	static LongTupleWord VERTEX_MASK;    // mask for bits in high word of tuple

	// the actual tuple
	LongTuple tuple;
	ssize_t count;

	ssize_t GetMult() {
		return 0;
	}

	static void SetTupleSize(int ts) {
		tupleSize = ts;                       // # nucleotides in a tuple
		numWords1 = (ts-1) / LT_NucsPerWord;  // index of high order word
		numWords = numWords1 + 1;             // number of words

		// Mask bits in high order word
		VERTEX_MASK = LT_Mask_Safe(ts % LT_NucsPerWord, 1) - LT_1;
	}
	
	void ZeroTuple() {
		for (int i = 0; i < LongTupleNumWords; i++) {
			tuple[i] = 0;
		}
	}

	IntegralTuple() {
		ZeroTuple();
		count = 0;
	}

	ssize_t Length() {
		return tupleSize;
	}

	ssize_t Valid() {
		// This may only store valid tuples.
		return 1;
	}

	ssize_t SetMult(ssize_t count) {
		return 0;
	}
	
	ssize_t ReadLine(std::istream &in, ssize_t minMult =0) {
		std::string stringTuple;
		if (!(in >> stringTuple)) {
			std::cout << "Error reading tuple." << std::endl;
			exit(1);
		}
		if (tupleSize == -1) {
			// determine the tuple size from the word that was read.
			SetTupleSize(stringTuple.size());
		}
		StringToTuple(stringTuple);
		std::string line;
		std::getline(in, line);
		// Parse the multiplicity
		std::stringstream linestrm(line);
		ssize_t mult = 0;
		linestrm >> mult;
		if (mult >= minMult) 
			return 1;
		else
			return 0;
	}

	// Convert string to a bitwise representation, two bits per nucleotide.
	// Store the tuple so that lowest two bits correspond to the first
	// nucleotide, and bits 2*k-1 & 2*k-2 (k is the tuple size)
	// correspond to the last nucleotide.
	// If s is a valid string (A,C,G,T only):
  //    set tuple to the bitwise representation
  //    return true
  // If s is not a valid string (e.g., "N"'s in sequence):
  //    tuple is set to garbage
  //    return false
	// TODO: Some calls to this don't check the return status.  Fix them.
	bool StringToTuple(const unsigned char *s) {
		LongTupleWord nuc;
		ssize_t p;
		ZeroTuple();


		for (p = 0; p < IntegralTuple::tupleSize; p++) {
			nuc = unmasked_nuc_index[s[p]];
			if (nuc >= 4)
				return 0;
			tuple[LT_WordNum(p)] |= LT_Mask_Mod(p, nuc);
		}
		return 1;
	}

	unsigned char GetNucleotide(ssize_t nucPos) {
		return (tuple[LT_WordNum(nucPos)] >> LT_BitNum(nucPos)) & LT_3;
	}

	// Faster routine for special case of 0th nucleotide
	unsigned char GetNucleotide0() {
		return tuple[0] & LT_3;
	}

	bool StringToTuple(const std::string &s) {
		return StringToTuple((unsigned char*) s.c_str());
	}

	// compare this.tuple and rhs.tuple
	// return
	//   -1: this.tuple < rhs.tuple
	//    0: this.tuple == rhs.tuple
	//    1: this.tuple > rhs.tuple

	int cmp(const IntegralTuple &rhs) const {
		int i;
		for (i = numWords1; i>=0; i--) {
			if (tuple[i] == rhs.tuple[i])
				continue;
			if (tuple[i] < rhs.tuple[i])
				return -1;
			return 1; // tuple > rhs.tuple
		}
		return 0; // tuple == rhs.tuple
	}

	static int cmp(const LongTuple &a, const LongTuple &b) {
		int i;
		for (i = LongTupleNumWords - 1; i>=0; i--) {
			if (a[i] == b[i])
				continue;
			if (a[i] < b[i])
				return -1;
			return 1; // a > b
		}
		return 0; // a == b
	}

	bool operator<(const IntegralTuple &rhs) const {
		return cmp(rhs) < 0;
	}

	bool operator>(const IntegralTuple &rhs) const {
		return cmp(rhs) > 0;
	}
	
	bool operator==(const IntegralTuple &rhs) const {
		return cmp(rhs) == 0;
	}
	
	bool operator!=(const IntegralTuple &rhs) const { 
		return cmp(rhs) != 0;
	}
	
	IntegralTuple& operator=(const IntegralTuple &rhs) {
		if (this != &rhs) {
			int i;
			for (i = 0; i < numWords; i++) {
				tuple[i] = rhs.tuple[i];
			}
		}
		return *this;
	}

#if 0
	IntegralTuple& operator=(const LongTuple &rhsTuple) {
		if (this != &rhs) {
			int i;
			for (i = 0; i < numWords; i++) {
				tuple[i] = rhsTuple[i];
			}
		}
		return *this;
	}
#endif

	void CopyTuple(const IntegralTuple &rhs) {
		if (this != &rhs) {
			int i;
			for (i = 0; i < numWords; i++) {
				tuple[i] = rhs.tuple[i];
			}
		}
	}

	IntegralTuple &operator=(const std::string &tupleString) {
		this->StringToTuple(tupleString);
		return *this;
	}


	// TODO: optimize.  E.g., could make lookup tables to do a byte at a time;
	// would need several lookup tables due to four possible shifts
	void MakeRC(IntegralTuple &rc) {

		rc.ZeroTuple();
		LongTupleWord word, nuc;
		for (int i = 0,
					 rcPos = tupleSize - 1;
				 i < tupleSize;
				 i++, rcPos--) {
			if ((i % LT_NucsPerWord) == 0) {
				// TODO: colorspace version w/o ~
				word = ~tuple[LT_WordNum(i)];
			}
			nuc = word & LT_3;
			word >>= 2;
			rc.tuple[LT_WordNum(rcPos)] |= LT_Mask_Mod(rcPos,nuc);
		}

	}

	void ToString() const;

	void ToString(std::string &tupleString) const {
		tupleString = "";
		tupleString.reserve(tupleSize);
		tupleString.resize(tupleSize);
		LongTupleWord word;
		for (int i = 0; i < tupleSize; i++) {
			if ((i % LT_NucsPerWord) == 0) {
				word = tuple[LT_WordNum(i)];
			}
			tupleString[i] = nuc_char[word & LT_3];
			word >>= 2;
		}
	}

	friend std::ostream &operator<<(std::ostream &out, const IntegralTuple &tup) {
		std::string tupleString;
		tup.ToString(tupleString);
		out << tupleString << " 0" << std::endl;
		return out;
	}

	friend std::istream &operator>>(std::istream &in, IntegralTuple &tup) {
		tup.ReadLine(in);
		return in;
	}
	
	ssize_t IncrementMult() {
		// No-op
		return 0;
	}

};

class CountedIntegralTuple : public IntegralTuple {
 public:

 CountedIntegralTuple() : IntegralTuple() {
		count = 1;
	}

	ssize_t IncrementMult() {
		if (count < MAX_COUNT)
			return ++count;
		else
			return (ssize_t) count;
		return 1;
	}

	ssize_t GetMult() {
		return (ssize_t) count;
		//		return 1;
	}

	ssize_t SetMult(ssize_t mult) {
		if (mult < MAX_COUNT)
			count = mult;
		else 
			count = MAX_COUNT;
		return count;
	}

	// Advance count by mult, but don't go over MAX_COUNT
	// mult must be positive
	ssize_t AdvanceMult(ssize_t mult) {
		ssize_t newcount = count + mult;
		if (newcount < count     // integer overflow in +
				|| newcount > MAX_COUNT) {
			count = MAX_COUNT;
		} else {
			count = newcount;
		}
		return count;
	}

	CountedIntegralTuple& operator=(const CountedIntegralTuple &rhs) {
		if (this != &rhs) {
			this->CopyTuple(rhs);
			//			this->tuple = rhs.tuple;
			this->count = rhs.count;
		}
		return *this;
	}
};

class EdgePosIntegralTuple : public IntegralTuple {
 public:
	ssize_t edge;
	ssize_t pos;
 EdgePosIntegralTuple() : IntegralTuple() {
	}
	EdgePosIntegralTuple& operator=(const EdgePosIntegralTuple &rhs) {
		if (this != &rhs) {
			this->CopyTuple(rhs);
			//			this->tuple = rhs.tuple;
			this->edge  = rhs.edge;
			this->pos   = rhs.pos;
		}
		return *this;
	}
};


template<typename T_Tuple>
ssize_t LookupBinaryTuple(T_Tuple *tupleList, ssize_t listLength, T_Tuple &value) {

	ssize_t begin = 0, end = listLength;
	ssize_t cur = (end + begin) / 2;
	while (begin < end and tupleList[cur] != value) {
		if (value < tupleList[cur]) {
			end = cur;
		}
		else {
			begin = cur + 1;
		}
		cur = (end + begin) / 2;
	}
	if (begin == end or tupleList[cur] != value)
		return -1;
	else
		return cur;

	/*

	//UNUSED// T_Tuple *ptr;
	ptr = std::lower_bound(tupleList, tupleList + listLength, value);
	//ptr = std::find(tupleList, tupleList + listLength, value);
	if (ptr != tupleList + listLength and *ptr == value)
		return ptr - tupleList;
	else
		return -1;
	*/
}


template<typename T_Tuple> 
ssize_t LookupBinaryTuple(std::vector<T_Tuple> &tupleList, T_Tuple &value) {
	typename std::vector<T_Tuple>::iterator tupleIt;
	tupleIt = std::lower_bound(tupleList.begin(), tupleList.end(), value);
	if (tupleIt != tupleList.end() and *tupleIt == value)
		return (tupleIt - tupleList.begin());
	else
		return -1;
}

template<typename T_Tuple>
void ReadBinaryTupleList(std::string &fileName, 
												 T_Tuple **tupleList,
												 ssize_t &listLength, ssize_t minMult = 0,
												 std::ostream &report = std::cout
												 ) {
	std::ifstream in;
	openck(fileName, in, std::ios::in | std::ios::binary, report);
	// read the number of tuples

	if (minMult == 0) {
		ssize_t tupleSize_SSZT;
		in.read((char*) &tupleSize_SSZT, sizeof(tupleSize_SSZT));
		T_Tuple::SetTupleSize(tupleSize_SSZT);

		in.read((char*) &listLength, sizeof(ssize_t));
		*tupleList = new T_Tuple[listLength];
		in.read((char*) *tupleList, sizeof(T_Tuple) * listLength);
	}
	else {
		ssize_t tupleSize_SSZT;
		in.read((char*) &tupleSize_SSZT, sizeof(tupleSize_SSZT));
		T_Tuple::SetTupleSize(tupleSize_SSZT);

		T_Tuple tuple;
		ssize_t totalTuples;
		listLength = 0;
		in.read((char*) &totalTuples, sizeof(ssize_t));
	 
		while(in) {
			if (in.read((char*) &tuple, sizeof(T_Tuple))) {
				if (tuple.GetMult() >= minMult) {
					//tupleList.push_back(tuple);
					listLength++;
				}
			}
		}
		if (listLength == 0 ) {
			*tupleList = NULL;
			return;
		}
		cout << "Placing " << listLength << " " << T_Tuple::tupleSize << "-mers into a dictionary." << endl;
		in.clear();
		in.seekg(std::ios::beg);
		*tupleList = new T_Tuple[listLength];
		ssize_t i = 0;

		in.read((char*) &tupleSize_SSZT, sizeof(tupleSize_SSZT));
		in.read((char*) &totalTuples, sizeof(ssize_t));
		while(in) {
			if (in.read((char*) &tuple, sizeof(CountedIntegralTuple))) {
				if (tuple.GetMult() >= minMult) {
					//tupleList.push_back(tuple);
					//				std::string tupleStr;
					//				tuple.ToString(tupleStr);
					//					(*tupleList)[i].tuple = tuple.tuple;
					(*tupleList)[i].CopyTuple(tuple);
					(*tupleList)[i].SetMult(tuple.GetMult());
					i++;
				}
			}
		}
	}
}


ssize_t ReadMultBoundedBinaryTupleList(std::string &fileName,
																	 ssize_t minMult,
																	 std::vector<CountedIntegralTuple> &tupleList);

template<typename T_Tuple>
void WriteBinaryTupleList(std::string &fileName,
													T_Tuple *tupleList,
													 ssize_t &listLength) {

	std::ofstream out;
	openck(fileName, out, std::ios::out | std::ios::binary);


	ssize_t tupleSize_SSZT = T_Tuple::tupleSize;
	out.write((const char*) &tupleSize_SSZT, sizeof(tupleSize_SSZT));
	out.write((const char*) &listLength, sizeof(ssize_t));
	out.write((const char*) tupleList, sizeof(T_Tuple) * listLength);
	out.close();
}


// The dictionary lookup table is an index based on the final IndexN
// nucleotides of a tuple, so it fits in a LongTuple.
// Since it is an index into a memory-based table, it is also a size_t.
// So it must fit in both, but it is really a size_t.

template<typename T_Tuple, ssize_t IndexN>
class DictionaryTupleList{
 public:
	_SSZT_ *startIndex, *endIndex;
	_SSZT_ indexLength;
	T_Tuple *list;
	_SSZT_ listLength;
	LongTupleWord mask;
	ssize_t offset;
	void Initialize(T_Tuple *listp, _SSZT_ listLengthp) {
		list = listp;
		listLength = listLengthp;
		this->CreateIndex();
	}

	size_t ComputeIndex(LongTuple tuple) {
		size_t index;
		// TODO: could be optimzed by defining as static class variables
		int lowWordNum = LT_WordNum(offset);
		int lowBitNum = LT_BitNum(offset);


		/*********************************************************
		 * Case: index is wholly within the most significant word
		 *********************************************************/

#if LongTupleNumWords == 1
		index = (tuple[T_Tuple::numWords1] >> lowBitNum);
		return index;
#else
		if (lowWordNum == T_Tuple::numWords1) {
			index = (tuple[T_Tuple::numWords1] >> lowBitNum);
			return index;
		}
#endif

		/*********************************************************
		 * Case: index is split over two or more words
		 *********************************************************/

		/* grab bits from top word */
		index = tuple[T_Tuple::numWords1];

		/* grab bits from middle words */

		// If sizeof(LongTupleWord) >= sizeof(size_t)
		// then lowWordNum = numWords1 or numWords1-1,
		// so the following loop logically could never execute.
		// It only has a chance to execute if LongTupleWord is smaller than size_t.
    // So don't compile it unless it's potentially needed.
#ifdef LongTupleWord_less_size_t
		for (int j = T_Tuple::numWords1 - 1; j > lowWordNum; j--) {
			index <<= LT_BitsPerWord;
			index |= (size_t) tuple[j];
		}
#endif

		/* grab bits from low word */
		size_t rem = (tuple[lowWordNum] >> lowBitNum);
		index <<= (LT_BitsPerWord - lowBitNum);
		index |= rem;
		return index;

	}

	void IndexList() {	Pause();
		if (IndexN == 0) {
			indexLength = 0;
			return;
		}
		assert(2*IndexN < SIZE_BITS);

		indexLength = 1;
		indexLength <<= (2*IndexN);
		startIndex = new ssize_t[indexLength];
		endIndex   = new ssize_t[indexLength];
		
		// Create a mask to grab the upper bits of the tuple.
		//
		offset = T_Tuple::tupleSize - IndexN;

		// TODO: new code doesn't use mask, but does have other constants that could
		// be precomputed.
		mask = ~(LT_Mask_Safe(offset,1) - LT_1);
		
		size_t upperIndex;
		
		_SSZT_ indexPos, listPos;
		indexPos = listPos = 0;
		
		// Loop invariant: 
		//  index[indexPos] is the first value i such that (list[i] & mask) == indexPos <= i

		_SSZT_ listStart = -1;
		listPos = indexPos = 0;
		while ((listPos < listLength) and 
					 (indexPos < indexLength)) {
			
			upperIndex = ComputeIndex(list[listPos].tuple);

			while (indexPos < indexLength and indexPos != upperIndex ) {
				// The index does not exist in the list
				assert(indexPos < indexLength);
				startIndex[indexPos] = -1;
				endIndex[indexPos]   = -1;
				indexPos++;
			}
			if (indexPos >= indexLength)
				break;

			// the beginning index where list[i] & mask == upperIndex
			listStart = listPos;

			while (indexPos == upperIndex) {
				// Advance forward in the list.
				listPos++;
				if (listPos >= listLength)
					break;
				upperIndex = ComputeIndex(list[listPos].tuple);
			}
			// The array from 'indexStart .. indexPos-1' has indexPos as it's uper 2*n bits.
			startIndex[indexPos] = listStart;
			endIndex[indexPos]   = listPos;
			
			// Searched all entries beginning with 'indexPos'
			indexPos++;
		}
	}
	
	ssize_t DictLookupBinaryTuple(T_Tuple& query) {

		size_t index = ComputeIndex(query.tuple);
		// Look if this tuple is indexed or not
		if (startIndex[index] == -1) {
			return -1;
		}

		ssize_t offsetIndex;
		offsetIndex = LookupBinaryTuple((T_Tuple*) &list[startIndex[index]], 
																		endIndex[index] - startIndex[index], query); 
		if (offsetIndex == -1)
			return -1;
		else
			return offsetIndex + startIndex[index];

	}
	
	void InitFromFile(std::string &fileName, ssize_t minMult = 0,
										std::ostream &report = std::cout
										) {
		ReadBinaryTupleList(fileName, &list, listLength, minMult, report);
		IndexList();
	}
};


typedef std::vector<CountedIntegralTuple> CountedIntegralTupleList;
typedef std::vector<EdgePosIntegralTuple> EdgePosIntegralTupleList;

typedef DictionaryTupleList<CountedIntegralTuple,9> CountedIntegralTupleDict;
typedef DictionaryTupleList<IntegralTuple,9> IntegralTupleDict;


template <typename T_Tuple>
void ForwardNuc(T_Tuple curVertex, char nextNuc, 
								T_Tuple &nextVertex) {
	assert(nextNuc < 4);

	int i;
	LongTupleWord nuc;
	nextVertex.ZeroTuple();
	for (i = 0; i < IntegralTuple::numWords1; i++) {
		nuc = curVertex.tuple[i+1] & LT_3;
		nuc = LT_Mask_Mod(LT_NucsPerWord-1, nuc);
		nextVertex.tuple[i] = nuc | (curVertex.tuple[i] >> 2);
	}
	nuc = LT_Mask_Mod(IntegralTuple::tupleSize-1, nextNuc);
	nextVertex.tuple[IntegralTuple::numWords1] =
		nuc | (curVertex.tuple[IntegralTuple::numWords1] >> 2);
}

template <typename T_Tuple>
void BackwardsNuc(T_Tuple &curVertex, char nextNuc,
									T_Tuple &prevVertex) {
	int i;
	int highBitNum = 2*(LT_NucsPerWord-1);
	LongTupleWord nuc = nextNuc;

	prevVertex.ZeroTuple();
	for (i = 0; i < IntegralTuple::numWords; i++) {
		prevVertex.tuple[i] = nuc | (curVertex.tuple[i] << 2);
		nuc = curVertex.tuple[i] >> highBitNum;
	}
	prevVertex.tuple[IntegralTuple::numWords1] &= IntegralTuple::VERTEX_MASK;
}


template <typename T_Tuple>
void MutateTuple(T_Tuple curTuple, unsigned char newNuc, ssize_t pos, T_Tuple &nextTuple) {
	int wordIndex = LT_WordNum(pos);
	int nucPos = pos % LT_NucsPerWord;

	curTuple.tuple[wordIndex] =
		curTuple.tuple[wordIndex] & ~LT_Mask(nucPos,3);

	nextTuple.CopyTuple(curTuple);
	nextTuple.tuple[wordIndex] =
		curTuple.tuple[wordIndex] | LT_Mask(nucPos,newNuc);
}



template <typename TupleType>
class UpdateIntegralTupleCountFunctor {
 public:
	ssize_t operator()(TupleType &tuple) {
		return tuple.IncrementMult();
	}
};




// TODO: Is IntegralTupleHashValueFunctor in use any more?
// The hash function seems to grab lower bits instead of upper bits.
template<ssize_t T_HashLength>
class IntegralTupleHashValueFunctor {
 public:
	LongTupleWord mask;

	IntegralTupleHashValueFunctor() {
		mask = LT_Mask_Safe(T_HashLength,1) - LT_1;
	}

	ssize_t operator()(const LongTuple longTuple) const {

		if (2*T_HashLength <= LT_BitsPerWord) {
			return (ssize_t) (longTuple[0] & mask);
		}
		size_t index = 0;
		for (int i = 0; i < LT_WordNum(T_HashLength-1); i++) {
			size_t subWord = longTuple[i];
			subWord <<= LT_BitsPerWord * i;
			index |= subWord;
		}
		index &= (size_t(1) << (2*T_HashLength)) - size_t(1);
		return (ssize_t) index;


	}

	// TODO: merge with above operator()
	ssize_t operator()(const IntegralTuple &tuple) const {
		if (2*T_HashLength <= LT_BitsPerWord) {
			return (ssize_t) (tuple.tuple[0] & mask);
		}

		// Use size_t instead of ssize_t so that 0's get shifted in,
		// 0's get padded in when doing type conversion from LongTupleWord
		// to size_t (in case LongTupleWord is smaller), etc.
		size_t index = 0;
		for (int i = 0; i <= LT_WordNum(T_HashLength-1); i++) {
			size_t subWord = tuple.tuple[i];
			subWord <<= LT_BitsPerWord * i;
			index |= subWord;
		}

		// Grab the bottom 2*T_HashLength bits
		// Note that the final word ored in above may have higher numbered bits.
		index &= (size_t(1) << (2*T_HashLength)) - size_t(1);
		return (ssize_t) index;
	}

	ssize_t MaxHashValue() {
		return ssize_t(1) << (2*(T_HashLength));
	}
};

																	 
#endif
