/***************************************************************************
 * Title:          MaximizeStrip.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef _MAXIMIZE_STRIP_
#define _MAXIMIZE_STRIP_

#include "StripGen.h"

void MaximizeStrips(T_Strips *strips, ssize_t window, ssize_t threshold);
void MinThreadingPath(T_Strips *&strips, ssize_t window);
#endif
