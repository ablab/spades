/***************************************************************************
 * Title:          modify_tuple.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "TupleDef.h" 
#include <iostream>

static ssize_t nuc_rc[] = {2,3,0,1};

Tuple snp(Tuple index, ssize_t n, ssize_t v, ssize_t k) {
  Tuple mask = 1, temp;
  mask <<= 2* (k - n - 1);

  temp = index >> 2 * (k - n);  // save the left hand side of the index
  temp <<= 2; // make room for new nt
  temp += v;  // append the new nt.
  temp = temp << 2 * (k - n - 1);   // make room for rest of sequence
  temp += index % mask;
  return temp;
}

Tuple insert(Tuple orig, ssize_t pos, ssize_t nuc, ssize_t tuple_len) {
  // parameters:
  // index:  the 4-bit representation of a sequence string to insert a
  // nucleotide.
  // n:      the position to insert.
  // v:      the value (nucleotide) to insert.
  // k:      the length of the index.

  Tuple mask = 1, temp;
  // must treat 1 as long before left shifting it.
  mask = 1;
  mask <<= (tuple_len*2);
  temp = (orig >> (pos*2));
  temp <<= 2;
  temp +=nuc;
  temp <<= (pos*2);
  temp += (orig % ( 1 << (pos*2)));
  if (!mask) {
    std::cout << "ins: mask is 0 " << tuple_len*2 << std::endl;
    exit(1);
  }
  temp = temp % mask;
  return temp;
}

Tuple del(Tuple orig, ssize_t pos) {
  Tuple temp;
  Tuple mask = 1;
  mask <<= (pos*2);
  temp = orig;
  temp >>= ((pos+1)*2);
  temp <<= (pos*2);
  if (!mask) {
    std::cout << "del: mask is : "<< (pos *2)<< std::endl;
    exit(1);
  }
     
  temp += orig % mask;
  return temp;
}

Tuple ReverseComplement(Tuple orig, ssize_t tuple_len) {
  ssize_t i;
  Tuple result = 0;
  ssize_t nuc;
  Tuple mask;
  mask = 1 << 2;
  for (i = 0; i < tuple_len; i++ ) {
    nuc = orig % mask;
    result <<= 2;
    result += nuc_rc[nuc];
    orig >>= 2;
  }
  return result;
}


ssize_t getNuc(Tuple index, ssize_t n, ssize_t k) {
  index >>= ((k- 1) - n) * 2;
  return index % 4;
}


Tuple subset(Tuple index, ssize_t lb, ssize_t ub, ssize_t k) {

  // mask off lower bound
  Tuple mask = 1;
  mask <<= 2*lb;
  index %= mask;

  // shift off upper bound
  index >>= 2*(k-ub);

  return index;
}

