#ifndef _PEAK_DETECTOR_PARAMETERS_HH_
#define _PEAK_DETECTOR_PARAMETERS_HH_
#include <iostream>

using namespace std;

class PeakDetectorParameters {
  double minCoverage, maxCoverage, binWidth;
  double snRatio;
public:
  PeakDetectorParameters( double minCov, double maxCov, double width, double ratio )
    : minCoverage(minCov), maxCoverage(maxCov), binWidth(width), snRatio(ratio){}
  bool checkParameters() const;
  double getMinCoverage() const { return minCoverage; }
  double getMaxCoverage() const { return maxCoverage; }
  double getBinWidth() const { return binWidth; }
  double getSnRatio() const { return snRatio; }
};

#endif // _PEAK_DETECTOR_PARAMETERS_HH_
