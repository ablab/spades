/***************************************************************************
 * Title:          BlockLib.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/22/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef _BLOCK_LIB_
#define _BLOCK_LIB_

#include <stdlib.h>
#include <sys/types.h>

void GetOffsetPosition(ssize_t rStart, ssize_t rEnd, ssize_t rPos,
		      ssize_t qStart, ssize_t qEnd, ssize_t &qPos);

#endif
