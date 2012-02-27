/***************************************************************************
 * Title:          compcommon.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef _COMP_COMMON
#define _COMP_COMMON
void ParseFileName(std::string fileName, 
		   std::string &refName, 
		   std::string &qryName,
		   std::string &seqName);

#endif
