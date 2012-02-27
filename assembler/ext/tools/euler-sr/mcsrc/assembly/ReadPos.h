/***************************************************************************
 * Title:          ReadPos.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/26/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef READ_POS_H_
#define READ_POS_H_

#include "SimpleSequence.h"
#include "assert.h"
#include <memory.h>
#include <string.h>

class ReadPos {
 public:
  ssize_t read;
  ssize_t pos;
  static SimpleSequenceList *sequences;
  static ssize_t hashLength;

  ReadPos() { read = -1; pos = -1;}
  ReadPos(ssize_t r, ssize_t p) { read = r; pos = p;}
  bool operator==(const ReadPos &comp) const {
    assert(read >= 0 and  read < (ssize_t) sequences->size());
    assert(comp.read >= 0 and comp.read < (ssize_t) sequences->size());
    assert(pos + hashLength <= (*sequences)[read].length);
    assert(comp.pos + hashLength <= (*sequences)[comp.read].length);
    return (memcmp((const char*) &(*sequences)[read].seq[pos], 
									 (const char*) &(*sequences)[comp.read].seq[comp.pos],
									 hashLength) == 0);
  }
  bool operator!=(const ReadPos &comp) const {
    return ! (*this == comp);
  }
  ReadPos& operator=(const ReadPos &val) {
		if (this != &val) {
			read = val.read;
			pos  = val.pos;
		}
		return *this;
  }
	bool operator<(const ReadPos &val) const {
		return strncmp((const char*)&(*sequences)[read].seq[pos],
									 (const char*)&(*sequences)[val.read].seq[val.pos], hashLength) < 0;
	}
	
  friend std::ostream &operator<<(std::ostream &out, const ReadPos &rp) {
    out << rp.read << " " << rp.pos;
    return out;
  }
  friend std::istream &operator>>(std::istream &in, ReadPos &rp) {
    in >> rp.read >> rp.pos;
    return in;
  }
};

class CompareReadPos {
 public:
	ssize_t operator()(const ReadPos &val, const char* seq) {
		return strncmp((const char*) &(*val.sequences)[val.read].seq[val.pos],
									 seq, val.hashLength) < 0;
	}
};

typedef std::vector<ReadPos> ReadPositions;


class CountedReadPos : public ReadPos {
public:
  ssize_t count;

  CountedReadPos() : ReadPos() {
    count = 1;
  }
  CountedReadPos & operator=(const CountedReadPos &copy) {
		if (this != &copy) {
			count = copy.count;
			read  = copy.read;
			pos   = copy.pos;
		}
    return *this;
  }
	friend std::istream &operator>>(std::istream &in, CountedReadPos &rp) {
		in >> (ReadPos&) rp >> rp.count;
		return in;
	}
	friend std::ostream &operator<<(std::ostream &out, const CountedReadPos &rp) {
		out << (ReadPos&) rp << " " << rp.count;
		return out;
	}
};

class UpdateFunctor {
public:
  ssize_t operator()(CountedReadPos &p) {
    p.count++;
		return p.count;
  }
};


template<typename L>
class CompareTuples {
 public:
  L *sequencesPtr;
  ssize_t length;
  ssize_t operator()(const ReadPos &a, const ReadPos &b) {
    return (strncmp((const char*) &((*sequencesPtr)[a.read].seq[a.pos]), 
										(const char*) &((*sequencesPtr)[b.read].seq[b.pos]),
										length) < 0);
  }
};

void PrintTuple(SimpleSequenceList &sequences, ReadPos &readPosition, 
								int tupleSize, std::ostream &out = std::cout);

/*
void PrintTuple(ssize_t tuple, int tupleSize, std::ostream &out);
*/

#endif
