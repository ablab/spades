/***************************************************************************
 * Title:          Alignment.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/26/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_



class Alignment {
public:
  ssize_t refPos;
  ssize_t qryPos;
  ssize_t refEnd, qryEnd;
  ssize_t *locations;
  ssize_t alignmentLength;
  Alignment() {
    refPos = qryPos = -1;
    refEnd = -1; qryEnd = -1;
    locations = NULL;
  }
  Alignment &operator=(const Alignment &copy) {
		if (this != &copy) {
			refPos = copy.refPos;
			qryPos = copy.qryPos;
			refEnd = copy.refEnd;
			refPos = copy.refPos;
			locations = copy.locations;
		}
		return *this;
  }
};
 
class T_Alignment : public Alignment {
public:
  ssize_t refEnum;
  ssize_t qryEnum;
  ssize_t length;
  ssize_t *enumerations;
  ssize_t level;
  ssize_t hashLength;
  T_Alignment *parent;
  T_Alignment *next;
  T_Alignment *child;

  T_Alignment() : Alignment() {
    refPos = 0;
    qryPos = 0;
    refEnum = -1;
    qryEnum = -1;
    length = 0;
    level  = 1;
    hashLength = 0;
    enumerations = NULL;
    locations    = NULL;
    next         = NULL;
    child        = NULL;
    parent       = NULL;
  }
};

void ComputeAlignmentStatistics(ssize_t *alignment, ssize_t alignLength, char *qry, char* ref,
																ssize_t &nMatch, ssize_t &nMismatch, ssize_t &nGapQry, ssize_t &nGapRef);

#endif
