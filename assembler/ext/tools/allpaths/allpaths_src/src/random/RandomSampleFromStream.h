///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef RANDOM_SAMPLE_FROM_STREAM
#define RANDOM_SAMPLE_FROM_STREAM

#include "Vec.h"

/// Take a random sample of a certain length from a stream of unknown size.
///
/// \class RandomSampleFromStream
///
/// For a sample of length N, we can show that all events in a stream of 
/// length m have the same probability of being taken if we do the following:
/// - We take sample k with probability N/k
/// - if we have taken sample k, we replace one of the N samples at random
/// with k.
///
/// This is a functor class : construct it, and then call operator() with
/// the current sample from the stream.
///
/// RandomSampleFromStreamN1: sample as RandomSampleFromStream, but specialized to
/// the case N = 1.

template<class T>
class RandomSampleFromStream {
private:
  vec<T> sample_;
  const int N_;
  longlong k_;

public:
  RandomSampleFromStream(int N=2000): sample_(N), N_(N), k_(0) {}

  void operator()(const T & item) {
    ++k_;
    if (k_ <= N_) sample_[k_-1] = item;
    else if (drand48() < double(N_)/k_) {
      int pos = int(floor(drand48() * N_));
      sample_[pos] = item;
    }
  }

  const vec<T> & Sample() { return sample_; }
};

template<class T>
class RandomSampleFromStreamN1 {
private:
  T sample_;
  unsigned int k_;

public:
  RandomSampleFromStreamN1( ): k_(0) {}

  bool Exists( ) { return k_ != 0; }

  void operator()(const T & item) {
    if (++k_ == 1) sample_ = item;
    else if (drand48() < 1.0/k_) sample_ = item;
  }

  const T& Sample() { return sample_; }
};

#endif // RANDOM_SAMPLE_FROM_STREAM
