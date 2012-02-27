/***************************************************************************
 * Title:          bitio.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
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
#ifndef BITIO_H
#define BITIO_H

#include <cstdio>
//#include <inttypes.h> // needs to be specially loaded in word.h
#include "word.h"

class bitio {
  FILE *f;
  word_t c;
  int  ph;
  bool writ;

 public:
  bool flush() {
    if (ph && writ) { 
      if (!fwrite(&c,sizeof(word_t),1,f)) return false;
    }
    ph = 0;
    c = 0;
    return true;
  }
  ~bitio() { flush(); }
  bitio(FILE *F,const char *p): f(F),c(0),ph(0),writ(false) {
    for(; *p; ++p) writ = writ || (*p == 'w') || (*p == 'a');
  };

  inline bool printbit(bool b) {
    c |= (word_t(b)<<ph);
    ph = (ph+1)%64;
    if (!ph) {
      if (!fwrite(&c,sizeof(word_t),1,f)) return false;
      c = 0;
    }
    return true;
  }
  inline bool printword(word_t w) {
    if (ph)
      c = (c << (64-ph)) | (w >> ph);
    else
      c = w;
    if (!fwrite(&c,sizeof(word_t),1,f))
      return false;

    c = w;
    return true;
  }

  inline bool scanbit(bool &b) {
    if (!ph && !fread(&c,sizeof(word_t),1,f)) return false;
    b = (c&( word_LIT(1) <<ph));
    ph = (ph+1)%64;
    return true;
  }
  inline bool scanword(word_t &w) {
    if (ph)
      w = (c << (64-ph));
    else
      w = 0;
    if (!fread(&c,sizeof(word_t),1,f)) return false;
    w |= (c >> ph);
    return true;
  }

  inline FILE *file() {return f;}
};
  
#endif
