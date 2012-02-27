/***************************************************************************
 * Title:          Net.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "Net.h"

ssize_t Net::fill = 0;
ssize_t Net::gap  = 1;

ssize_t Net::top = 0;
ssize_t Net::syn = 1;
ssize_t Net::inv = 2;
ssize_t Net::nonSyn = 3;


Net::Net() {
  level = -1;
  chainClass = -1;
  tStart = qStart = -1;
  tSize = qSize = -1;
  chrom = "";
  orientation =-1;
  id = -1;
  score = 0;
  ali = -1;
  qFar = -1;
  qOver = -1;
  qDup = -1;
  type= -1;
  tN = qN = tR = qR = tNewR = qNewR = tOldR = qOldR = tTrf = qTrf = -1;
}
