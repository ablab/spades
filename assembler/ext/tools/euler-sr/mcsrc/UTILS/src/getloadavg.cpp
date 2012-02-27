/***************************************************************************
 * Title:          getloadavg.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
int main() {
  double loads[3];
  ssize_t result;
  result = getloadavg(loads, 3);
  if (result < 0) {
    printf("error getting load\n");
    exit(1);
  }
  printf("%.2g %.2g %.2g\n", loads[0], loads[1], loads[2]);
  return 0;
}
