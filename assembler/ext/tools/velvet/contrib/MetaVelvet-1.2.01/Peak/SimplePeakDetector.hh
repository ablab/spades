#ifndef _SIMPLE_PEAK_DETECTOR_HH_
#define _SIMPLE_PEAK_DETECTOR_HH_
#include <math.h>
#include <vector>
#include <map>
#include "PeakDetectorParameters.hh"
#include "PeakDetector.hh"

using namespace std;

class SimplePeakDetector : public PeakDetector {
  map<double, double> data;
  double minCoverage, maxCoverage, binWidth;
  double snRatio;
  void initializeCounts();
  void setCounts( MetaGraph* graph );
  double getNoiseCutoff() const;
  double searchFirstValley() const;
  void setCount( double binVal, double count ){ data.insert( map<double, double>::value_type(binVal, count) ); }
  double getCount( double binVal ) const { return data.find(binVal)->second; }
  void addCount( double binVal, double count ){ data.find(binVal)->second += count; }
public:
  SimplePeakDetector( const PeakDetectorParameters* param ) 
    : minCoverage(param->getMinCoverage()), maxCoverage(param->getMaxCoverage()), binWidth(param->getBinWidth()), 
      snRatio(param->getSnRatio()){
    initializeCounts();
  }
  int detectCoveragePeaks( MetaGraph* graph, double* coveragePeaks, double* coverageBoundaries );
};

#endif // _META_HISTO_HH_
