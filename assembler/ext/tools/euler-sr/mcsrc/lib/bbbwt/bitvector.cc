/*
** Copyright 2004 Ross Lippert
** 
** This file is part of bbbwt.
** 
** bbbwt is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 2 of the License, or
** (at your option) any later version.
** 
** bbbwt is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with bbbwt; if not, write to the Free Software
** Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
**
*/
#include "bitvector.h"

void bitvector::fdump(FILE *f) const {
  for(unsigned i=0; i < size(); i++) {
    fputc( (char)('0'+(*this)[i]) , f);
  }
}

void bitvector::fprint(FILE *f) const {
  unsigned char c;
  
  for(unsigned i=7; i < size(); i+=8) {
    c = 0;
    c = c | ((*this)[i-0]);
    c = c | ((*this)[i-1]<<1);
    c = c | ((*this)[i-2]<<2);
    c = c | ((*this)[i-3]<<3);
    c = c | ((*this)[i-4]<<4);
    c = c | ((*this)[i-5]<<5);
    c = c | ((*this)[i-6]<<6);
    c = c | ((*this)[i-7]<<7);
    fputc( c , f);
  }
}

bool bitvector::fscan(FILE *f) {
  unsigned char c;
  int ci=0;
  for(unsigned i=7; i < size() && (ci = fgetc(f))!=EOF ; i+=8) {
    c = ci;
    (*this)[i-0] = (c)    & ((unsigned char)1);
    (*this)[i-1] = (c>>1) & ((unsigned char)1);
    (*this)[i-2] = (c>>2) & ((unsigned char)1);
    (*this)[i-3] = (c>>3) & ((unsigned char)1);
    (*this)[i-4] = (c>>4) & ((unsigned char)1);
    (*this)[i-5] = (c>>5) & ((unsigned char)1);
    (*this)[i-6] = (c>>6) & ((unsigned char)1);
    (*this)[i-7] = (c>>7) & ((unsigned char)1);
  }
  return ci != EOF;
}

