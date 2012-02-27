/***************************************************************************
 * Title:          MCHsp.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef MCHsp_H_
#define MCHsp_H_

#include <vector>

class MCHsp {
public:
  ssize_t refStart, refEnd, qryStart, qryEnd;
  double score;
  ssize_t *alignment;
  MCHsp(const MCHsp &copy) {
    refStart = copy.refStart;
    refEnd   = copy.refEnd;
    qryStart = copy.qryStart;
    qryEnd   = copy.qryEnd;
    score = copy.score;
  }
  MCHsp(ssize_t refStartP, ssize_t refEndP, 
	ssize_t qryStartP, ssize_t qryEndP) {
    refStart = refStartP;
    refEnd   = refEndP;
    qryStart = qryStartP;
    qryEnd   = qryEndP;
  }

  MCHsp() {
    refStart = refEnd = qryStart  = qryEnd = -1;
    score = 0.0;
  }
  MCHsp &operator=(const MCHsp &hsp) {
		if (this != &hsp) {
			refStart = hsp.refStart;
			refEnd   = hsp.refEnd;
			qryStart = hsp.qryStart;
			qryEnd   = hsp.qryEnd;
			score    = hsp.score;
		}
		return *this;
  }
};



class MCHspIntervalList {
public:
  std::vector<MCHsp*> hsps;
  
  MCHsp* Find(MCHsp *hsp, double ovpRatio);
  MCHsp* Insert(MCHsp *hsp, double ovpRatio);
};

/*
  typedef RBTreeNode<MCHsp> HSPTreeNode;
  typedef RBTree<MCHsp> HSPMap;
*/
#endif
