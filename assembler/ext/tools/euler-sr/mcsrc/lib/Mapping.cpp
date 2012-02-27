/***************************************************************************
 * Title:          Mapping.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "Mapping.h"


void PrintGlue(ssize_t *locations, ssize_t length, std::ostream& out) {
  ssize_t i;
  for (i = 0; i < length; i++) {
    if (locations[i] != -1)
      out << i << " @;" << locations[i] << " !;" << std::endl;
  }
}


void PrintAlignmentCompressed(ssize_t *locations, ssize_t length, std::ostream &out) {
  //UNUSED+// ssize_t j;
  ssize_t i ;
  i = 0;
  out << length << std::endl;
  //UNUSED+// ssize_t endRef,endQry;
  ssize_t startRef, startQry  ;
  while (i < length) {
    // skip past unaligned sequences
    while (locations[i] == -1)
      i++;
    startRef = i;
    startQry = locations[i];
    while (i < length - 1  && locations[i] == locations[i+1]-1) 
      i++;

    out << startRef << "\t" << startQry << "\t" << i - startRef + 1;
  }
}

void ReadGlue(ssize_t *&locations, ssize_t &length, std::istream& in) {
  // Count the number of glue
  


}
