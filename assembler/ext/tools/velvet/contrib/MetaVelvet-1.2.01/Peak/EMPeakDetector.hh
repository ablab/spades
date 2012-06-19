#ifndef _EM_PEAK_DETECTOR_HH_
#define _EM_PEAK_DETECTOR_HH_
#include <math.h>

#include "PeakDetector.hh"

#define EM_PEAK_DETECTOR_MAX_NUM_CLASS PEAK_DETECTOR_MAX_NUM_COVERAGE_PEAKS
#define EM_PEAK_DETECTOR_MAX_HISTO_LENGTH 10000

class EMPeakDetector : public PeakDetector {
  int fineboundary( double* histo, int nowboundary );
  int myEM( double *group, long double *newave, long double *newweight, int count, int groupNum, int classNum, long double prob[][EM_PEAK_DETECTOR_MAX_NUM_CLASS] );
  int judge_neighbor( int small, int large, int *index, int *classNum, long double *tmpave,long double *tmpweight, long double *ave, long double *weight, double n_parameter, int tmpindex );
  void detectPeakPands( double *smoothhisto, int xMax, int xMin, int binWidth, double thresHeight, double *listPeakPands);
  void smoothingHisto( double *newhisto, double *smoothhisto, int xMax, int xMin, int binWidth, int widthMovAve );
  void weightedHisto( double *histo, double *newhisto, int xMax, int xMin, int binWidth );
  void setWidthByxMax( int xMax, int *listWidth );
  int setXMax( int xMax, int binWidth );
public:
  int detectCoveragePeaks( MetaGraph* graph, double coveragePeaks[EM_PEAK_DETECTOR_MAX_NUM_CLASS], double coverageBoundaries[EM_PEAK_DETECTOR_MAX_NUM_CLASS] );
};

#endif // _EM_PEAK_DETECTOR_HH_
