/***************************************************************************
 * Title:          Inversion.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef _INVERSION_H_
#define _INVERSION_H_

#include <string>
#include <map>
#include <vector>



class Inversion {
public:
  std::string sequence;
  std::string species;
  ssize_t tStart, tEnd, qStart, qEnd;
};

typedef std::vector<Inversion*> InversionVector;
typedef std::map<std::string, InversionVector > SpeciesInversions;
typedef std::map<std::string, SpeciesInversions > InversionMap;


#endif
