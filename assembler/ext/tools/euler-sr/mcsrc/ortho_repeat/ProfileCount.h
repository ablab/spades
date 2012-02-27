/***************************************************************************
 * Title:          ProfileCount.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef PROFILE_COUNT_H_
#define PROFILE_COUNT_H_

#include "mctypes.h"



class ProfileCount {
public:
  ssize_t nChars, length;
  IntMatrix counts;
  IntVector total;
  ProfileCount(ssize_t _nchars, ssize_t _length) : nChars(_nchars), length(_length) {
    profile.resize(nChars);
    ssize_t i, j;
    for (i = 0; i < nChars; i++) {
      counts[i].resize(length);
      for (j =0 ; j < length; j++) 
	counts[i][j] = 0;
    }
    total.resize(length);
    for (j =0 ; j < length; j++) 
      total[j] = 0;
  }
  void Set(ssize_t pos, ssize_t value);
};  


#endif 
