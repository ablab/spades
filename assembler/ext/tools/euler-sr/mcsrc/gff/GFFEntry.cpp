/***************************************************************************
 * Title:          GFFEntry.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "GFFEntry.h"


bool GFFEntry::ParsePALSLine(std::istream &in) {
  if (!in.good() or in.eof())
    return false;

  std::string junk;
  // sample input line
  //  chr3  pals  hit  10863665  10864065  310  +  .  Target  chr1  17255 17657;  maxe  0.058
  if (!(in >> refSeq >> source >> feature >> refStart >> refEnd 
	>> score >> strand >> frame >> junk 
	>> qrySeq >> qryStart >> qryEnd 
	>> junk >> junk >> errorRatio)) return false;
  return true;
}
