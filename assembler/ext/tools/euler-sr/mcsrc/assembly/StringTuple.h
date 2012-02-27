/***************************************************************************
 * Title:          StringTuple.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  02/27/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef STRING_TUPLE_H_
#define STRING_TUPLE_H_

#include "DNASequence.h"
#include "Tuple.h"
#include "SeqUtils.h"

class StringTuple {
 public:
	static int tupleSize;
	char *s;
	ssize_t GetMult() { 
		return 0;
	}

	StringTuple() {
		s = NULL;
	}

	StringTuple(char* sp) {
		s = sp;
	}

	ssize_t Length() {
		return tupleSize;
	}
	char * ToString() {
		return (char*) s;
	}
	
	void assign(char *str) {
		s = str;
	}

	ssize_t IncrementMult() {
		return 0;
	}

	static void SetTupleSize(int ts) {
		tupleSize = ts;
	}
		
	ssize_t ReadLine(std::ifstream &in, ssize_t minMult=0) {
		//UNUSED// ssize_t mult;
		if (tupleSize == -1) {
			std::string temp;
			in >> temp;
			tupleSize = temp.size();
			if (s != NULL) delete[] s;
			s = new char[tupleSize];
			strncpy(s, temp.c_str(), tupleSize);
		}
		else {
			if (s == NULL) 
				s = new char[tupleSize];
			in >> s;
			std::string line;
			// pitch the remainder of the line since it may have a multiplicity
			// or other junk.
			std::getline(in, line); 
			std::stringstream linestrm(line);
			ssize_t mult = 0;
			linestrm >> mult;
			if (mult >= minMult) 
				return 1;
			else
				return 0;
		}
	}

	ssize_t Valid() {
		ssize_t i;
		if (s == NULL)
			return 0;
		for (i = 0; i < tupleSize; i++) {
			if (numeric_nuc_index[(unsigned char) s[i]] >= 4)
				return 0;
		}
		return 1;
	}

	ssize_t GetHashValue(size_t &hashValue) {
		ssize_t i;
		if (s == NULL) {
			return -1;
		}
		hashValue = 0;
		// method 1
		for (i = 0; i < tupleSize; i++) {
			hashValue <<=2;
			if (unmasked_nuc_index[(unsigned char) s[i]] >= 4)
				return -1;
			else
				hashValue += unmasked_nuc_index[(unsigned char) s[i]];
			/*			
			if (i < tupleSize-1)
			hashValue <<=2;*/
		}

		/*
		size_t h2 = 0;
		for (i = 0; i < tupleSize; i++) {
			h2 <<=2;
			if (unmasked_nuc_index[(unsigned char) s[i]] >= 4)
				return -1;
			else
				h2 += unmasked_nuc_index[(unsigned char) s[i]];
		}
		
		if (h2 != hashValue) {
			std::cout << "error: " << h2 << " != " << hashValue << std::endl;
			}*/
		return 1;
	}

	friend std::istream& operator>>(std::istream &in, StringTuple &tuple) {
		std::string templine;
		if (tupleSize == -1) {
			in >> templine;
			tupleSize = templine.size();
			if (tuple.s == NULL) {
				tuple.s = new char[tupleSize];
			}
			strncpy(tuple.s, templine.c_str(), tupleSize);
		}
		else {
			in >> tuple.s;
		}
		std::getline(in, templine);
		return in;
	}
	
	friend std::ostream &operator<<(std::ostream &out, const StringTuple &tuple) {
		if (tupleSize != -1 and tuple.s != NULL) {
			std::string strtup(tuple.s,tupleSize);
			out << strtup << std::endl;
		}
		return out;
	}
	
	bool operator<(const StringTuple &rhs) const {
		return strncmp(s, rhs.s, tupleSize) < 0;
	}
	
	bool operator>(const StringTuple &rhs) const {
		return strncmp(s, rhs.s, tupleSize) > 0;
	}
	bool operator==(const StringTuple &rhs) const {
		return strncmp(s, rhs.s, tupleSize) == 0;
	}
	bool operator!=(const StringTuple &rhs) const {
		return ! (*this==rhs);
	}
	unsigned char operator[](ssize_t pos) {
		assert(s != NULL);
		return s[pos];
	}
	StringTuple& copy(const StringTuple &rhs) {
		if (s == NULL) {
			s = new char[tupleSize];
		}
		memcpy(s, rhs.s, tupleSize);
		return *this;
	}

	StringTuple& operator=(const StringTuple &rhs) {
		if (this != &rhs) {
			this->copy(rhs);
		}
		return *this;
	}
			
};

//int  GetHashValue(DNASequence &seq, int pos, int length, StringTuple &tuple);


#endif
