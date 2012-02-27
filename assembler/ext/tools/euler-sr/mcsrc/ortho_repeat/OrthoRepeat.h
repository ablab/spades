/***************************************************************************
 * Title:          OrthoRepeat.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef ORTHO_REPEAT_H_
#define ORTHO_REPEAT_H_

#include <string>

class OrthoRepeat {
public:
  ssize_t refStart, refEnd;
  ssize_t refStrand;
  ssize_t qryStart, qryEnd;
  ssize_t qryStrand;
  std::string refType;
  std::string qryType;
};


#endif
