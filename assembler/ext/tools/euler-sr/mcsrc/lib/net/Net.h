/***************************************************************************
 * Title:          Net.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef NET_H_
#define NET_H_

#include <string>

class Net {
public:
  ssize_t level;
  ssize_t chainClass;
  ssize_t tStart, qStart;
  ssize_t tSize, qSize;
  std::string chrom;
  ssize_t orientation;
  ssize_t id;
  ssize_t score;
  ssize_t ali;
  ssize_t qFar;
  ssize_t qOver;
  ssize_t qDup;
  ssize_t type;
  ssize_t tN, qN, tR, qR, tNewR, qNewR, tOldR, qOldR, tTrf, qTrf;
  static ssize_t top, syn, inv, nonSyn;
  static ssize_t fill, gap;
  Net();
};

//myeostatin, prevent muscle growth

#endif
