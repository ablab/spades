/***************************************************************************
 * Title:          ValidatedInversion.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "ValidatedInversion.h"

void InversionMatrix::PrintMiserly(std::ostream &out){
  ssize_t l;
  for (l = 0; l < size(); l++ ) {
    loci[l]->PrintOrientations(out);
    out << std::endl;
  }
}
void InversionMatrix::PrintConcise(std::ostream &out) {
  ssize_t l;
  for (l = 0; l < size(); l++ ) {
    loci[l]->Print(out);
  }
}  

void InversionMatrix::Print(std::ostream &out) {
  if (title != "" )
    out << ">" << title << std::endl;
  PrintConcise(out);
}

void InversionMatrix::PrintOrdered(std::ostream &out, std::vector<std::string> &order) {
  ssize_t l;
  ssize_t s;

  for (s = 0; s < order.size(); s++) {
    for (l = 0; l < size(); l++ ) {
      if (loci[l]->species == order[s]) 
	loci[l]->Print(out);
    }
  }
}

void CopyInversionList(InversionList &src, InversionList &dest) {
  ssize_t i;
  dest.clear();
  for (i = 0; i < src.size(); i++ ) {
    dest.push_back(new InversionMatrix(*src[i]));
  }
}
