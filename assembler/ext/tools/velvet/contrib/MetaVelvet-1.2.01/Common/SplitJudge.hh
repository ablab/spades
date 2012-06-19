#ifndef _SPLIT_JUDGE_HH_
#define _SPLIT_JUDGE_HH_
#include <math.h>
#include "../VelvetAPI/VelvetGraph.hh"
#include "SplitStats.hh"

#define META_GRAPH_MAX_NUM_COVERAGE_PEAKS 100

using namespace std;

class SplitJudge {
  double repeatCoverageSD;
  size_t minSplitLength;
  bool flagUseConnections;
  size_t numValidConnections, numNoiseConnections;
  int numCoveragePeaks;
  double expectedCoverages[META_GRAPH_MAX_NUM_COVERAGE_PEAKS];
  SplitStatsWriter* statsWriter;
  bool flagReportDetail;
  bool isLocalRepeatStructure( Node* node ) const;
  bool isPassLengthFilter( Node* node ) const;
  bool isPassLengthFilter( const vector<Node*>& nodes ) const;
  bool isRepeatCoverageCondition( Node* node ) const;
  bool isConsistentConnections( Node* node, SplitStats* stats ) const;
  size_t calcNumConsistentConnections( Node* node ) const;
  size_t calcNumConsistentConnectionsForPrimaryNodes( Node* node ) const;
  size_t calcNumConsistentConnectionsForNonPrimaryNodes( Node* node ) const;
  size_t calcNumInConsistentConnections( Node* node ) const;
  size_t calcNumConnections( Node* node1, Node* node2 ) const;
  double getNearestPeak( double coverage ) const ;
public:
  SplitJudge( double rcSD, size_t len, bool flagUC, size_t numVC, size_t numNC ) 
    : repeatCoverageSD(rcSD), minSplitLength(len), 
      flagUseConnections(flagUC),numValidConnections(numVC),  numNoiseConnections(numNC), 
      statsWriter( new SplitStatsWriter() ) {}
  void setStatsWriter( const string& statsFileName, const string& detailFileName, bool flagReportDetail );
  void setExpectedCoverages( int num, const double* covs );
  void showExpectedCoverages() const;
  int getNumCoveragePeaks() const { return numCoveragePeaks; }
  double getExpectedCoverage( size_t i ) const { return expectedCoverages[i]; }
  bool getFlagUseConnections() const { return flagUseConnections; }
  bool isSplit( Node* node ) const;
};

#endif // _SPLIT_JUDGE_HH_
