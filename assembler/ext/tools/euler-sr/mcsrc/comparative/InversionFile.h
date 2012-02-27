/***************************************************************************
 * Title:          InversionFile.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef _INVERSION_FILE_H_
#define _INVERSION_FILE_H_

#include <iostream>
#include <fstream>

#include <string>

#include "Inversion.h"
#include "mctypes.h"

ssize_t ReadInversionFile(std::string &invFileName, 
		      InversionMap &invMap,
		      StringSet &species);


#endif
