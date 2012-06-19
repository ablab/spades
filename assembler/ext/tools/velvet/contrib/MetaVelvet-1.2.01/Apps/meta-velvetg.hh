#ifndef _META_VELVET_G_HH_
#define _META_VELVET_G_HH_
#include <stdio.h>
#include <string.h>
#include <string>
#include <iostream>
#include "../VelvetAPI/VelvetHeaders.hh"
#include "../VelvetAPI/VelvetUtils.hh"
#include "../Common/FileUtils.hh"
#include "../Common/MetaGraph.hh"
#include "../Peak/PeakDetectorParameters.hh"
#include "../Peak/PeakDetectorFactory.hh"
#include "FileNameUtils.hh"
#include "AppsDefs.hh"

#define META_VELVET_G_PARAMETER_VALUE_NOT_SPECIFIED -1
#define META_VELVET_G_SD_RATIO_FOR_GUESS 0.10
#define META_VELVET_G_COVERAGE_DELIM_CHAR "_"
#define META_VELVET_G_MAX_NUM_COVERAGE_PEAKS PEAK_DETECTOR_MAX_NUM_COVERAGE_PEAKS
#define META_VELVET_G_DEFAULT_MIN_PEAK_COVERAGE 0.0
#define META_VELVET_G_DEFAULT_MAX_PEAK_COVERAGE 500.0
#define META_VELVET_G_DEFAULT_HISTO_BIN_WIDTH 1.0
#define META_VELVET_G_DEFAULT_HISTO_SN_RATIO 10.0
#define META_VELEVTG_DEFAULT_MAX_CHIMERA_RATE 0.0
#define META_VELVETG_DEFAULT_REPEAT_COVERAGE_SD 0.10
#define META_VELVETG_DEFAULT_MIN_SPLIT_LENGTH 0
#define META_VELVETG_DEFAULT_NUM_VALID_CONNECTIONS 1
#define META_VELVETG_DEFAULT_NUM_NOISE_CONNECTIONS 0

using namespace std;

class MetaVelvetG {
  string prefix;
  bool flagEstimateCoverageCutoff;
  bool flagEstimateExpectedCoverage;
  bool flagEstimateExpectedCoverages;
  bool flagReadTracking;
  bool flagScaffolding;
  bool flagExportUnusedReads;
  bool flagExportFilteredNodes;
  bool flagExportAmos;
  bool flagExportAlignments;
  boolean flagMatePair[CATEGORIES];
  int accelerationBits;
  long minContigLength;
  long minContigKmerLength;
  double coverageMask;
  double coverageCutoff;
  double maxCoverageCutoff;
  double longCoverageCutoff;
  double expectedCoverage;
  double insertLength[CATEGORIES];
  double insertLengthSD[CATEGORIES];
  double longInsertLength;
  double longInsertLengthSD;
  bool flagDiscardChimera;
  int numCoveragePeaks;
  double maxChimeraRate;
  double repeatCoverageSD;
  double expectedCoverages[META_VELVET_G_MAX_NUM_COVERAGE_PEAKS];
  double coverageBoundaries[META_VELVET_G_MAX_NUM_COVERAGE_PEAKS];
  size_t minSplitLength;
  double minPeakCoverage, maxPeakCoverage, histoBinWidth, histoSnRatio;
  bool flagUseConnections;
  size_t numValidConnections, numNoiseConnections;
  bool flagReportSplitDetail;
  bool flagReportSubgraph;
  void setDefaultParameters();
  bool checkParameters() const;
  bool checkMetaHistoParameters() const;
  bool checkParameters( int argc, char* argv[] ) const;
  bool checkArgc( int argc ) const;
  bool checkPrefix( int argc, char* argv[] ) const;
  bool checkMandatoryFiles( string prefix ) const;
  bool checkCategory( int cat ) const;
  bool checkCoverageOrder() const;
  void printUsage() const;
  bool setParameters( int argc, char* argv[] );
  bool setParameter( char* name, char* val );
  bool setCoverageCutoff( char* val );
  bool setExpectedCoverage( char* val );
  bool setLongCoverageCutoff( char* val );
  bool setInsertLength( int cat, char* val );
  bool setInsertLengthSD( int cat, char* val );
  bool setLongInsertLength( char* val );
  bool setLongInsertLengthSD( char* val );
  bool setFlagReadTracking( char* val );
  bool setFlagScaffolding( char* val );
  bool setFlagExportFilteredNodes( char* val );
  bool setFlagExportAmos( char* val );
  bool setFlagExportAlignments( char* val );
  bool setMinContigLength( char* val );
  bool setCoverageMask( char* val );
  bool setAccelerationBits( char* val );
  bool setMaxBranchLength( char* val );
  bool _setMaxDivergence( char* val );
  bool setMaxGapCount( char* val );
  bool setMinPairCount( char* val );
  bool setMaxCoverage( char* val );
  bool setLongMultiCutoff( char* val );
  bool _setPairedExpFraction( char* val );
  bool setFlagExportUnusedReads( char* val );
  bool setFlagMatePair( int cat, char* val );
  bool setFlagDiscardChimera( char* val );
  bool setMaxChimeraRate( char* val );
  bool setExpectedCoverages( char* val );
  bool setRepeatCoverageSD( char* val );
  bool setMinSplitLength( char* val );
  bool setMinPeakCoverage( char* val );
  bool setMaxPeakCoverage( char* val );
  bool setHistoBinWidth( char* val );
  bool setHistoSnRatio( char* val );
  bool setFlagUseConnections( char* val );
  bool setNumNoiseConnections( char* val );
  bool setNumValidConnections( char* val );
  bool setFlagReportSplitDetail( char* val );
  bool setFlagReportSubgraph( char* val );
  string getInitialGraphFileName() const;
  void guessNonMandatoryParameters();
  void guessInsertLengthSD();
  static bool isGuessInsertLengthSD( double ave, double sd );
  static double guessInsertLengthSD( double ave );
  void guessCoverageFlags();
public:
  MetaVelvetG( int argc, char* argv[] );
  void constructPregraph() const;
  MetaGraph* loadGraph() const;
  void setInsertLengths( MetaGraph* metaGraph ) const;
  void setSplitJudge( MetaGraph* metaGraph ) const;
  void setSubgraphOptions( MetaGraph* metaGraph ) const;
  void estimateCoverage( MetaGraph* metaGraph );
  void estimateExpectedCoverages( MetaGraph* metaGraph );
  PeakDetectorParameters* getPeakDetectorParameters() const;
  void estimateExpectedCoverage( MetaGraph* metaGraph );
  void estimateCoverageCutoff( MetaGraph* metaGraph );
  void removeNodes( MetaGraph* metaGraph );
  void scaffolding( MetaGraph* metaGraph );
  void finalize( MetaGraph* metaGraph );
  void output( MetaGraph* metaGraph );
};

#endif // _META_VELVET_G_HH_
