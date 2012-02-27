/***************************************************************************
 * Title:          Profile.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef PROFILE_H_
#define PROFILE_H_

// Class for actg, - profile

#include <vector>

#include "mctypes.h"

class Profile {
public:
  ssize_t nChars, length;
  FloatMatrix profile;
  Profile(ssize_t _nchars, ssize_t _length) : nChars(_nchars), length(_length) {
    profile.resize(nChars);
    ssize_t i, j;
    for (i = 0; i < nChars; i++) {
      profile[i].resize(length);
      for (j =0 ; j < length; j++) 
	profile[i][j] = 0;
    }
  }
};  

#endif
