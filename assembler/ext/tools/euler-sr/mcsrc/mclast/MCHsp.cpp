/***************************************************************************
 * Title:          MCHsp.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include "MCHsp.h"

MCHsp* MCHspIntervalList::Find(MCHsp *hsp, double ovpRatio) {
  ssize_t i;

  /*
    std::cout << "looking for " << hsp->refStart << " " << hsp->refEnd 
    << " " << hsp->qryStart << " " << hsp->qryEnd << std::endl;
  */


  for (i = 0; i < hsps.size(); i++) {
    // find overlapping hsp
    /*    std::cout << "checking " << hsps[i]->refStart << " " << hsps[i]->refEnd 
	      << " " << hsps[i]->qryStart << " " << hsps[i]->qryEnd << std::endl;
    */
    if (((hsp->refStart >= hsps[i]->refStart and 
	  hsp->refStart <= hsps[i]->refEnd) or
	 (hsp->refEnd >= hsps[i]->refStart and
	  hsp->refEnd <= hsps[i]->refEnd) or 
	 (hsp->refStart <= hsps[i]->refStart and
	  hsp->refEnd >= hsps[i]->refEnd)) and 
	((hsp->qryStart >= hsps[i]->qryStart and 
	  hsp->qryStart <= hsps[i]->qryEnd) or
	 (hsp->qryEnd >= hsps[i]->qryStart and
	  hsp->qryEnd <= hsps[i]->qryEnd) or 
	 (hsp->qryStart <= hsps[i]->qryStart and
	  hsp->qryEnd >= hsps[i]->qryEnd))) {
      // Found both reference and query overlap.
      ssize_t refOvp, qryOvp;
      refOvp = std::min(hsp->refEnd, hsps[i]->refEnd) - 
	std::max(hsp->refStart, hsps[i]->refStart); 
      qryOvp = std::min(hsp->qryEnd, hsps[i]->qryEnd) - 
	std::max(hsp->qryStart, hsps[i]->qryStart); 
      if (hsp->refEnd == hsp->refStart) {
	std::cout << "error, hsp has 0 size " << std::endl;
	exit(1);
      }
      if (hsp->qryEnd == hsp->refStart) {
	std::cout << "error, hsp has 0 size " << std::endl;
	exit(1);
      }

      if (double(refOvp)/ double(hsp->refEnd - hsp->refStart) <= ovpRatio and
					double(qryOvp)/ double(hsp->qryEnd - hsp->qryStart) <= ovpRatio) {
	/*
	  std::cout << "found it " << std::endl;
	*/
	return hsps[i];
      }
    }
  }
  return NULL;
}

MCHsp* MCHspIntervalList::Insert(MCHsp *hsp, double ovpRatio) {
  MCHsp *existing;
  existing = Find(hsp, ovpRatio); 
  if ( existing == NULL )
    hsps.push_back(hsp);

  return existing;
}
