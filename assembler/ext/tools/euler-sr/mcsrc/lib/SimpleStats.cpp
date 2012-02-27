/***************************************************************************
 * Title:          SimpleStats.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "SimpleStats.h"
#include "compatibility.h"

void SeedRandom() {
	time_t t;
  srandom((unsigned) time(&t));
}

_LONG_ Random(_LONG_ randMax) {
  return (_LONG_) (randMax * (double(random()) / RANDOM_MAX));
}

double Uniform() {
  return double(random()) / RANDOM_MAX;
}
	  
double NormalRandom(double mean, double var) {

  //UNUSED+// double y2;
  double x1, x2, w, y1 ;
  // Do Box-Muller transformation to generate N(0,1)
  do {
    x1 = 2.0 * Uniform() - 1.0;
    x2 = 2.0 * Uniform() - 1.0;
    w = x1 * x1 + x2 * x2;
  } while ( w >= 1.0 );
  w = sqrt( (-2.0 * log( w ) ) / w );
  y1 = x1 * w;
  return y1*var + mean;
}
