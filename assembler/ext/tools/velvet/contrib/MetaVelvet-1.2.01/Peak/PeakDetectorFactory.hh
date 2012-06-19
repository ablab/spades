#ifndef _PEAK_DETECTOR_FACTORY_HH_
#define _PEAK_DETECTOR_FACTORY_HH_
#include "PeakDetector.hh"
#include "PeakDetectorParameters.hh"
#include "EMPeakDetector.hh"
#include "SimplePeakDetector.hh"

using namespace std;

class PeakDetectorFactory {
public:
  static PeakDetector* instantiatePeakDetector( const string& algo, const PeakDetectorParameters* params );
};

#endif // _PEAK_DETECTOR_FACTORY_HH_
