/***************************************************************************
 * Title:          ReadPosTuple.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/26/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifnedf READ_POS_TUPLE_H_
#define READ_POS_TUPLE_H_

#include "ReadPos.h"
#include "Tuple.h"

class ReadPosTuple : public CountedReadPos, public Tuple {
 public:
	static SimpleSequenceList *sequences;
	ssize_t mult;
	ssize_t GetMult() {
		return count;
	}
	char * GetString() {
		return (*sequences)[read].seq[pos];
bj	}
	void assign(char *str) { 
		std::cout << "cannot assign readpos tuple to a string " << std::endl;
		exit(1);
	}
	
	ssize_t IncrementMult() {
		count++;
	}
		
	static ssize_t ReadLine(std::ifstream &in, TupleStr &tuple) {};
	ssize_t Valid()= 0;
	bool operator<(const TupleStr &rhs) const =0;
	bool operator>(const TupleStr &rhs) const =0;
	bool operator==(const TupleStr &rhs) const =0;
	bool operator !=(const TupleStr &rhs) const =0;
	Tuple &operator=(const Tuple &rhs)=0;



}

#endif
