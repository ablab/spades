#ifndef _VELVET_GRAPH_HH_
#define _VELVET_GRAPH_HH_
#include <iostream>
#include "VUtils.hh"

using namespace std;

class VelvetGraph {
protected:
  ReadSet* sequences;
  Graph* graph;
  Connection** scaffold;
  boolean* dubious;
  long getMinimumContigKmerLength( long minContigLength ) const;
  void removeLowCoverageNodes( double coverageCutoff, double longCoverageCutoff, long minContigKmerLength, bool flagExportFilteredNodes, 
			       const string& lowCoverageContigsFileName );
  void _removeHighCoverageNodes( double maxCoverageCutoff, long minContigKmerLength, bool flagExportFilteredNodes, const string& highCoverageContigsFileName );
public:
  VelvetGraph( const string& seqFileName, const string& roadmapFileName, const string& pregraphFileName, bool flagReadTracking, int accelerationBits );
  VelvetGraph( const string& seqFileName, const string& graphFileName );
  Graph* getGraph(){ return graph; }
  void setInsertLength( int cat, double ave, double sd );
  long getNumNodes() const;
  Node* getNode( long nodeID ) const;
  double estimateExpectedCoverage( const string& prefix );
  double estimateCoverageCutoff( const string& prefix );
  void correct();
  void removeNodes( double coverageCutoff, double longCoverageCutoff, double maxCoverageCutoff, 
		    long minCnotigLength, bool flagExportFilteredNodes, 
		    const string& lowCoverageContigsFileName, const string& highCoverageContigsFileName );
  void scaffoldingWithSinglePeakMode( double expectedCoverage, boolean* flagMatePair, bool flagScaffolding );
  void finalize( double coverageCutoff, double longCoverageCutoff );
  void save( const string& graphFileName ) const;
  void save( const string& graphFileName, long minContigLength, double coverageMask ) const;
  void saveStats( const string& statsFileName ) const;
  void saveAlignments( const string& alignmentFileName, long minContigLength, const string& seqFileName ) const;
  void saveAmos( const string& amosFileName, long minContigLength ) const;
  void saveUnusedReads( const string& prefix, long minContigLength ) const;
};

#endif // _VELVET_GRAPH_HH_
