// Copyright (c) 2004 Broad Institute of MIT and Harvard
// 
// NormalRandom.h

#ifndef NORMAL_RANDOM_H
#define NORMAL_RANDOM_H

#include "random/Random.h"
#include "random/NormalDistribution.h"


// class NormalRandom
//
// Provides a source of random numbers with a normal distribution.
//

class NormalRandom
{
public:
  // standard constructor gives a standard normal distribution
  NormalRandom()
    : mean_(0.0), stddev_(1.0)
  {}

  // constructor for normal distribution with specified mean and standard deviation
  NormalRandom(double mean, double stddev)
    : mean_(mean), stddev_(stddev)
  {}

  // method for getting the next value of the normal random variable  
  double value() const;

private:
  double mean_;
  double stddev_;
};

// FastNormal initially makes a list of 10^6 evaluations of NormalRandom r( 0, 1 ),
// and for each call returns the next member of the list, cycling back to the
// beginning upon reaching the end.

double FastNormal( );

#endif
