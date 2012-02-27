/***************************************************************************
 * Title:          BlockLib.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "BlockLib.h"

// This file will probably be deprecated pretty soon by 
// functions written for BlockGraph, and functions inside block.
// 

#include <math.h>
#include <stdlib.h>

#include "compatibility.h"

void GetOffsetPosition(ssize_t rStart, ssize_t rEnd, ssize_t rPos,
		      ssize_t qStart, ssize_t qEnd, ssize_t &qPos) {

  ssize_t rLen = szabs(rEnd - rStart);
  
  ssize_t qLen = szabs(qEnd - qStart);
  
  qPos = (ssize_t) floor((double(rPos - rStart)/ rLen)* qLen) + qStart;
}
