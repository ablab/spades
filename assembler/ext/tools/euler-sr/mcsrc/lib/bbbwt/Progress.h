/***************************************************************************
 * Title:          Progress.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
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
#ifndef PROGRESS_H
#define PROGRESS_H

#include <cstdio>

template<class Num> class Progress {
  int pct;
  Num count;
 public:
  Num maximum;
  FILE *file;

  Progress(): pct(0), count(0), maximum(0), file(NULL) {}

  ~Progress() {}

  inline void reset() { pct = 0; count = 0;}

  inline void add(Num a) {
    count += a;
    if (file && (pct*maximum < 100*count)) report();
  }

  inline void set(Num a) {
    count = a;
    if (file && (pct*maximum < 100*count)) report();
  }

 private:
  void report() {
    while (pct*maximum < 100*count) {
      fprintf(file," %d%% complete\r",++pct);
      fflush(file);
    }
  }
};
  

#endif
