/***************************************************************************
 * Title:          StringMultTuple.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/26/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef STRING_MULT_TUPLE_H_
#define STRING_MULT_TUPLE_H_

#include "StringTuple.h"

class StringMultTuple : public StringTuple {
 public:
	ssize_t mult;
	ssize_t GetMult() const {
		return mult;
	}
	char *GetString() const {
		return s;
	}

 StringMultTuple() : StringTuple() {
		mult = 0;
	}

 StringMultTuple(char *sp) : StringTuple(sp){
		mult = 0;
	}
	
	StringMultTuple &operator=(const StringMultTuple &rhs) {
		if (this == &rhs)
			return *this;

		if (s == NULL)
			s = new char[tupleSize+1];
		strncpy(s, rhs.s, tupleSize);
		mult = rhs.GetMult();
		return *this;
	}
	
	void ToString(std::string &str) {
		str.copy(s, tupleSize);
	}

	ssize_t IncrementMult(ssize_t inc = 1) {
		return mult+= inc;
	}
		
	ssize_t ReadLine(std::ifstream &in, ssize_t minMult = 0) {	
		if (tupleSize == -1) {
			std::string tuple;
			in >> tuple;
			tupleSize = tuple.size();
			if (s == NULL) {
				s = new char[tupleSize+1];
				s[tupleSize] = 0;
			}
			strncpy(s, tuple.c_str(), tupleSize);
			in >> mult;
		}
		else {
			if (s == NULL) {
				s = new char[tupleSize+1];
				s[tupleSize] = 0;
			}
			in >> s >> mult;
		}
		std::string remainder;
		std::getline(in, remainder);
		if (mult >= minMult)
			return 1;
		else
			return 0;
	}

	bool operator<(const StringMultTuple &rhs) const {
		return (strncmp(s, rhs.s, tupleSize) < 0);
	}

	bool operator>(const StringMultTuple &rhs) const {
		return (strncmp(s, rhs.s, tupleSize) > 0);
	}
	
	bool operator==(const StringMultTuple &rhs) const {
		return (strncmp(s, rhs.s, tupleSize) == 0);
	}
	
	bool operator !=(const StringMultTuple &rhs) const {
		return (!(*this == rhs));
	}

	friend std::ostream &operator<<(std::ostream &out, const StringMultTuple &rhs) {
		out << rhs.s << " " << rhs.mult << std::endl;
		return out;
	}
	friend std::istream &operator>>(std::istream &in, StringMultTuple &rhs) {
		if (rhs.s == NULL)
			rhs.s = new char[StringMultTuple::tupleSize];
		in >> rhs.s >> rhs.mult;
		return in;
	}
	void CleanUp() {
		if (s != NULL ){
			delete[] s;
		}
	}
	char operator[](ssize_t pos) {
		return s[pos];
	}
};


#endif
