/***************************************************************************
 * Title:          BlastHSP.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/28/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef BLAST_HSP_H_
#define BLAST_HSP_H_

#include <string>
#include "Alignment.h"

class BlastHSP : public Alignment {
 public:
  std::string query;
	std::string sbjct;
  double score;
  long double eValue;
  ssize_t strand;
	double identity;
	ssize_t length;
  BlastHSP() : Alignment() {
    score = 0; eValue = 0; strand = 0;
    identity = 0;
		query = "";
		sbjct = "";
  }

	// TODO: Fix sparc compiler warning
	// "warning: base class `class Alignment' should be explicitly initialized in the copy constructor"
  BlastHSP(const BlastHSP &hsp) {
		*this = hsp;
  }

  BlastHSP& operator=(const BlastHSP &hsp) {
		if (this != &hsp) {
			sbjct = hsp.sbjct;
			query = hsp.query;
			score = hsp.score;
			eValue = hsp.eValue;
			strand = hsp.strand;
			identity = hsp.identity;
			(Alignment(*this)) = (Alignment) hsp;
			alignmentLength = hsp.alignmentLength;
			refEnd = hsp.refEnd;
			qryEnd = hsp.qryEnd;
			refPos = hsp.refPos;
			qryPos = hsp.qryPos;
			locations = hsp.locations;
		}
		return *this;
  }
};



#endif
