/***************************************************************************
 * Title:          GapFunctions.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <string>

ssize_t MuchGreater(ssize_t a, ssize_t b, double ratio);
ssize_t GetSpeciesNameFromFile(std::string fileName, std::string &speciesName);
ssize_t GetRegionNameFromFile(std::string fileName, std::string &regionName);
