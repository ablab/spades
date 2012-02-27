/***************************************************************************
 * Title:          SimpleStats.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef _STATS_H_
#define _STATS_H_

#include "compatibility.h"

#include <vector>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>

void SeedRandom();

_LONG_ Random(_LONG_ randMax);
double Uniform();
	  
double NormalRandom(double mean, double var);

template <typename T>
T Max(std::vector<T> &values) {
  ssize_t i;
	if (values.size() == 0) {
		return -99999999;
	}
  T max = values[0];
  for (i = 1; i < values.size(); i++) {
    if (max < values[i])
			max = values[i];
	}
  
  return max;
}

template <typename T>
T Min(std::vector<T> &values) {
  ssize_t i;
	if (values.size() == 0) {
		return 99999999;
	}

  T min = values[0];
  for (i = 1; i < values.size(); i++) {
    if (min > values[i])
			min = values[i];
	}
  return min;
}

template <typename T>
double Integrate(std::vector<T> bins, double &width) {
  ssize_t i;
  ssize_t nbins;
  double area;
  nbins = bins.size();
  area = 0;
  for (i = 0; i < nbins; i++)
    area += bins[i] * width;
  return area;
}

template <typename T>
void Bin(std::vector<T> &values, std::vector<ssize_t> &bins, ssize_t nbins, T min, T max) {
  ssize_t i;
  bins.resize(nbins);

  // initialize bins
  for (i = 0; i < nbins; i++) {
    bins[i] = 0;
  }

  // calculate bins
  ssize_t binIndex;
  double ratio;
  //UNUSED// double lenth;
  double length;
  length = double(max) - double(min);
  for (i = 0; i < values.size(); i++) {
    if (values[i] == max)
      bins[nbins-1]++;
    else {
      ratio = (values[i] - min) / length;
      binIndex = (ssize_t) floor(ratio*nbins);
      bins[binIndex]++;
    }
  }
}


template <typename T>
void GetTail(std::vector<T> &values, std::vector<ssize_t> &tailIndices, double tailPct) {
  ssize_t i;
  double width, area;
  T min, max;
  std::vector<ssize_t> bins;
  ssize_t nbins = 100;

  min = Min(values);
  max = Max(values);
  Bin(values, bins, 100, min, max);
  width = (double(max) - double(min)) / nbins;
  area = Integrate(bins, width);
  double tailArea = 0.0;
  i = nbins  - 1;
  // how much of the tail makes tailPct volume? 
  while ( i > 0 && (tailArea / area) < tailPct) {
    tailArea += width*bins[i];
    i--;
  }
  double minTailValue;
  
  minTailValue = i*width + min;
  for (i = 0; i < values.size(); i++) 
    if (values[i] > minTailValue) 
      tailIndices.push_back(i);
}


template <typename T>
void GetMeanVar(std::vector<T> &values, T &mean, T &var) {

  ssize_t i;
	// TODO: should we use double or T for computing
	// (as opposed to returning) total, totalSSE, mean, var ?
  T total;
  T totalSSE;

  if (values.size() == 0) {
    mean = 0;
    var = -1; // n'est pas possible!
		return;
  }

  total = 0;
  totalSSE = 0;
  ssize_t length;
  for (i = 0; i < values.size(); i++) {
    total += values[i];
    totalSSE += values[i] * values[i];
  }

  mean = total / values.size();
  var  =
		T(values.size()) / (values.size()-1)
		* (totalSSE / values.size() - mean*mean);
}

#endif
