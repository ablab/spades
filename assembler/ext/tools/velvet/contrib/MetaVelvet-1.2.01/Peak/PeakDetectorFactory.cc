#include "PeakDetectorFactory.hh"

PeakDetector* PeakDetectorFactory::instantiatePeakDetector( const string& algo, const PeakDetectorParameters* params ){
  if( algo == "simple" ){
    return new SimplePeakDetector( params );
  } else if( algo == "em" ){
    return new EMPeakDetector();
  }
  return new SimplePeakDetector( params );
}
