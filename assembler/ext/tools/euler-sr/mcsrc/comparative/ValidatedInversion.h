/***************************************************************************
 * Title:          ValidatedInversion.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef VALIDATED_INVERSION_H_
#define VALIDATED_INVERSION_H_

#include <string>
#include <map>
#include <vector>
#include <set>
#include <iostream>
#include <iomanip>
#include <assert.h>
#include "mctypes.h"



class ValidatedInversion {
public:
  BitVector orient;
  ssize_t number;
  std::string species;
  ssize_t size() {
    return orient.size();
  }
  void PrintOrientations(std::ostream &out=std::cout) {
    ssize_t i;
    for (i = 0; i < orient.size(); i++) {
      out << std::setw(2) << orient[i];
    }
  }
  void Print(std::ostream &out=std::cout) {
    out << species << "\t";
    PrintOrientations(out);
    out << std::endl;
  }
  ValidatedInversion() { }
  ValidatedInversion(ValidatedInversion *copy) {
    assert(copy != NULL);
    *this = *copy;
  }

  ValidatedInversion &operator=(const ValidatedInversion& copy) {
		if (this != &copy) {
			species = copy.species;
			number  = copy.number;
			orient.resize(copy.orient.size());
			ssize_t i;
			for (i=  0; i < orient.size(); i++ ) {
				orient[i] = copy.orient[i];
			}
		}
  }
  ssize_t operator[](ssize_t i) {
    return orient[i];
  }
  void Erase(ssize_t pos) {
    assert(pos <orient.size());
    orient.erase(orient.begin() + pos);
  }
};

class InversionMatrix {
public:
  std::string title;
  std::vector<ValidatedInversion*> loci;
  ssize_t size() { return loci.size(); }
  ssize_t numSpecies() {
    if (loci.size() == 0)
      return 0;
    else 
      return loci[0]->orient.size();
  }
  void Print(std::ostream &out);
  void PrintMiserly(std::ostream &out);
  void PrintConcise(std::ostream &out);
  void PrintOrdered(std::ostream &out, std::vector<std::string> &order);
  InversionMatrix() { }
  InversionMatrix(InversionMatrix &copy) {
    ssize_t i;
    title = copy.title;
    for (i = 0; i < copy.size(); i++ ){
      loci.push_back(new ValidatedInversion(copy.loci[i]));
    }
  }
  void Erase(ssize_t index) {
    assert(index < size());
    loci.erase(loci.begin() + index);
  }
  const ValidatedInversion* & operator[](const ssize_t &index) {
    assert(index < size());
    return loci[index];
  }
};


typedef std::vector<InversionMatrix*> InversionList;
typedef std::map<ssize_t, ssize_t> BinIndexMap;
typedef std::vector<std::vector<double> > SimMat;

void CopyInversionList(InversionList &src, InversionList &dest);
#endif
