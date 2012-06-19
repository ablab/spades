#ifndef _PEAK_DETECTOR_HH_
#define _PEAK_DETECTOR_HH_
#include "../Common/MetaGraph.hh"

#define PEAK_DETECTOR_MAX_NUM_COVERAGE_PEAKS 100

using namespace std;

class PeakDetector {
public:
  virtual int detectCoveragePeaks( MetaGraph* graph, double coveragePeaks[PEAK_DETECTOR_MAX_NUM_COVERAGE_PEAKS], 
				   double coverageBoundaries[PEAK_DETECTOR_MAX_NUM_COVERAGE_PEAKS] ) = 0;
};

#endif // _PEAK_DETECTOR_HH_
