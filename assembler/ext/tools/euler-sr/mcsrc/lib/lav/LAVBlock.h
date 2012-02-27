/***************************************************************************
 * Title:          LAVBlock.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/26/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef LAVBLOCK
#define LAVBLOCK

#include <iostream>
#include <ostream>
#include <vector>
class LAVBlock {
public:
  ssize_t score;
  ssize_t refBegin, qryBegin;
  ssize_t refEnd, qryEnd;
  ssize_t enumeration;
  ssize_t strand;
  static ssize_t ref;
  static ssize_t qry;
  LAVBlock(ssize_t rb, ssize_t re, ssize_t qb, ssize_t qe, ssize_t st, ssize_t sc) {
    refBegin = rb;
    refEnd   = re;
    qryBegin = qb;
    qryEnd   = qe;
    score    = sc;
    strand   = st;
  }

  LAVBlock(ssize_t rb, ssize_t re, ssize_t qb, ssize_t qe, ssize_t st, ssize_t sc, ssize_t id) {
    refBegin = rb;
    refEnd   = re;
    qryBegin = qb;
    qryEnd   = qe;
    score    = sc;
    strand   = st;
    refALBegin.push_back(rb);
    refALEnd.push_back(re);
    qryALBegin.push_back(qb);
    qryALEnd.push_back(qe);
    alIdentity.push_back(id);
  }
  void AddInterval(ssize_t rb, ssize_t re, ssize_t qb, ssize_t qe, ssize_t id) {
    refALBegin.push_back(rb);
    refALEnd.push_back(re);
    qryALBegin.push_back(qb);
    qryALEnd.push_back(qe);
    alIdentity.push_back(id);
  }
  void CopyIntervals(ssize_t start, LAVBlock &block, ssize_t intvStart, ssize_t intvEnd) {
    ssize_t i;
    for (i = intvStart; i < intvEnd; i++ ) {
      refALBegin[i-intvStart + start] = block.refALBegin[i];
      refALEnd[i-intvStart + start]   = block.refALEnd[i];
      qryALBegin[i-intvStart + start] = block.qryALBegin[i];
      qryALEnd[i-intvStart+ start]   = block.qryALEnd[i];
    }
  }

  

  LAVBlock() {
    score = 0;
    refBegin = refEnd = qryBegin = qryEnd = -1;
    strand = -1;
    enumeration = 0;
  }
  ssize_t IntvBegin(ssize_t i, ssize_t sequence ) {
    if (sequence == 0) 
      return refALBegin[i];
    else
      return qryALBegin[i];
  }
  ssize_t IntvEnd(ssize_t i, ssize_t sequence ) {
    if (sequence == 0)
      return refALEnd[i];
    else
      return qryALEnd[i];
  }

  ssize_t begin(ssize_t sequence) {
    if (sequence == 0) 
      return refBegin;
    else
      return qryBegin;
  }
  
  ssize_t end(ssize_t sequence) {
    if (sequence == 0)
      return refEnd;
    else
      return qryEnd;
  }

  bool operator<(LAVBlock &comp) const {
    return qryBegin < comp.qryBegin;
  }
  void Resize(ssize_t size) {
    refALBegin.resize(size);
    refALEnd.resize(size);
    qryALBegin.resize(size);
    qryALEnd.resize(size);
    alIdentity.resize(size);
  }

  ssize_t size() { return refALBegin.size(); }
  std::vector<ssize_t> refALBegin, refALEnd, qryALBegin, qryALEnd, alIdentity;
  std::ostream &TabbedPrintBlock(std::ostream &out, ssize_t &chainId, ssize_t strand, 
				 ssize_t refSeqId = 0, ssize_t qrySeqId=0) {
    ssize_t i;
    for (i = 0; i < refALBegin.size(); i++) {
      out << 0 << "\t" 
	  << refALBegin[i] << "\t" << refALEnd[i] << "\t" 
	  << qryALBegin[i] << "\t" << qryALEnd[i] << "\t" << strand << "\t"
	  << chainId << "\t" << refSeqId << "\t" << qrySeqId << std::endl;
      ++chainId;
    }
    return out;
  }

  std::ostream &PrintBlock(std::ostream &out) {
    out << "a {" << std::endl;
    out << "  s " << score << std::endl;
    out << "  b " << refBegin << " " << qryBegin << std::endl;
    out << "  e " << refEnd << " " << qryEnd << std::endl;
    ssize_t i;
    for (i = 0; i < refALBegin.size(); i++) {
      out << "  l " << refALBegin[i] << " " << qryALBegin[i] << " " << refALEnd[i] 
	  << " " << qryALEnd[i]   << " " << alIdentity[i] << std::endl; 
    }
    out << "}" << std::endl;
    return out;
  }
};

#endif
