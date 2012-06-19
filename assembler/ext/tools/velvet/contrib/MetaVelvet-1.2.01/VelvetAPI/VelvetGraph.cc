#include "VelvetGraph.hh"

VelvetGraph::VelvetGraph( const string& seqFileName, const string& roadmapFileName, const string& pregraphFileName, bool flagReadTracking, int accelerationBits ){
  sequences = importReadSet( VUtils::getCharFileName(seqFileName) );
  convertSequences( sequences );
  graph = importPreGraph( VUtils::getCharFileName(pregraphFileName), sequences, VUtils::getCharFileName(roadmapFileName), (bool)flagReadTracking, accelerationBits );
}

VelvetGraph::VelvetGraph( const string& seqFileName, const string& graphFileName ){
  sequences = importReadSet( VUtils::getCharFileName(seqFileName) );
  convertSequences( sequences );
  graph = importGraph( VUtils::getCharFileName(graphFileName) );
}

void VelvetGraph::setInsertLength( int cat, double ave, double sd ){
  setInsertLengths( graph, (Category)cat, ave, sd );
}

long VelvetGraph::getNumNodes() const {
  return (long)nodeCount(graph);
}

Node* VelvetGraph::getNode( long nodeID ) const {
  return getNodeInGraph( graph, nodeID );
}

double VelvetGraph::estimateExpectedCoverage( const string& prefix ){
  return estimated_cov( graph, VUtils::getCharFileName(prefix) );
}

double VelvetGraph::estimateCoverageCutoff( const string& prefix ){
  return estimated_cov( graph, VUtils::getCharFileName(prefix) ) / 2;
}

void VelvetGraph::correct(){
  cout << "[VelvetGraph] " << " === Correct graph === " << endl;
  ShortLength* sequenceLengths = getSequenceLengths( sequences, getWordLength(graph) );
  correctGraph( graph, sequenceLengths, sequences->categories );
}

void VelvetGraph::removeNodes( double coverageCutoff, double longCoverageCutoff, double maxCoverageCutoff, 
			       long minContigLength, bool flagExportFilteredNodes, 
			       const string& lowCoverageContigsFileName, 
			       const string& highCoverageContigsFileName ){
  long minContigKmerLength = getMinimumContigKmerLength( minContigLength );

  cout << "[VelvetGraph] " << " === Remove low coverage nodes ===" << endl;
  removeLowCoverageNodes( coverageCutoff, longCoverageCutoff, minContigKmerLength, flagExportFilteredNodes, lowCoverageContigsFileName );

  cout << "[VelvetGraph] " << " === Remove high coverage nodes ===" << endl;
  _removeHighCoverageNodes( maxCoverageCutoff, minContigKmerLength, flagExportFilteredNodes, highCoverageContigsFileName );

  cout << "[VelvetGraph] " << " === Clip tips hardly ===" << endl;
  clipTipsHard(graph);
  
  if (sequences->readCount > 0 && sequences->categories[0] == REFERENCE){
    cout << "[VelvetGraph] " << " === Removing low coverage arcs ===" << endl;
    removeLowArcs(graph, coverageCutoff);
  }
}

void VelvetGraph::removeLowCoverageNodes( double coverageCutoff, double longCoverageCutoff, long minContigKmerLength, bool flagExportFilteredNodes, 
					  const string& lowCoverageContigsFileName ){
  dubious = removeLowCoverageNodesAndDenounceDubiousReads( graph, coverageCutoff, sequences, flagExportFilteredNodes, 
							   (Coordinate)minContigKmerLength, VUtils::getCharFileName(lowCoverageContigsFileName) );
  removeLowLongCoverageNodesAndDenounceDubiousReads( graph, longCoverageCutoff, sequences, dubious, flagExportFilteredNodes, 
						     (Coordinate)minContigKmerLength, VUtils::getCharFileName(lowCoverageContigsFileName) );
}

void VelvetGraph::_removeHighCoverageNodes( double maxCoverageCutoff, long minContigKmerLength, bool flagExportFilteredNodes, 
					    const string& highCoverageContigsFileName ){
  removeHighCoverageNodes( graph, maxCoverageCutoff, (Coordinate)minContigKmerLength, flagExportFilteredNodes, 
			   VUtils::getCharFileName(highCoverageContigsFileName) );
}

long VelvetGraph::getMinimumContigKmerLength( long minContigLength ) const {
  if( minContigLength < 2 * getWordLength(graph) ){
    return getWordLength(graph);
  }
  return minContigLength - getWordLength(graph) + 1;
}

void VelvetGraph::scaffoldingWithSinglePeakMode( double expectedCoverage, boolean* flagMatePair, bool flagScaffolding ){
  if(expectedCoverage <= 0){ return; }

  cout << "[VelvetGraph] " << " === Rock Bank  ===" << endl;
  readCoherentGraph( graph, isUniqueSolexa, expectedCoverage, sequences );

  cout << "[VelvetGraph] " << " === Create read paring array ===" << endl;
  createReadPairingArray(sequences);

  cout << "[VelvetGraph] " << " === Detach dubious reads ===" << endl;
  detachDubiousReads( sequences, dubious );

  cout << "[VelvetGraph] " << " === Activate gap markers ===" << endl;
  activateGapMarkers(graph);

  cout << "[VelvetGraph] " << " === Scaffolding ===" << endl;
  for ( int pebbleRounds=pairedCategories(sequences)+1 ; pebbleRounds>0 ; --pebbleRounds){
    exploitShortReadPairs( graph, sequences, dubious, flagMatePair, (boolean)flagScaffolding );
  }
}

void VelvetGraph::finalize( double coverageCutoff, double longCoverageCutoff ){
  if( dubious ){ free(dubious); }
  concatenateGraph( graph );
  removeLowCoverageReferenceNodes( graph, coverageCutoff, longCoverageCutoff, sequences );
}

void VelvetGraph::save( const string& graphFileName ) const {
  exportGraph( VUtils::getCharFileName(graphFileName), graph, sequences->tSequences );
}

void VelvetGraph::save( const string& graphFileName, long minContigLength, double coverageMask ) const {
  long minContigKmerLength = getMinimumContigKmerLength( minContigLength );
  ShortLength* sequenceLengths = getSequenceLengths( sequences, getWordLength(graph) );
  exportLongNodeSequences( VUtils::getCharFileName(graphFileName), graph, minContigKmerLength, sequences, sequenceLengths, coverageMask); 
}

void VelvetGraph::saveStats( const string& statsFileName ) const {
  displayGeneralStatistics( graph, VUtils::getCharFileName(statsFileName), sequences );
}

void VelvetGraph::saveAlignments( const string& alignmentFileName, long minContigLength, const string& seqFileName ) const {
  long minContigKmerLength = getMinimumContigKmerLength( minContigLength );
  exportLongNodeMappings( VUtils::getCharFileName(alignmentFileName), graph, sequences, minContigKmerLength, VUtils::getCharFileName(seqFileName) );
}

void VelvetGraph::saveAmos( const string& amosFileName, long minContigLength ) const {
  long minContigKmerLength = getMinimumContigKmerLength( minContigLength );
  exportAMOSContigs( VUtils::getCharFileName(amosFileName), graph, minContigKmerLength, sequences );
}

void VelvetGraph::saveUnusedReads( const string& prefix, long minContigLength ) const {
  long minContigKmerLength = getMinimumContigKmerLength( minContigLength );
  exportUnusedReads( graph, sequences, minContigKmerLength, VUtils::getCharFileName(prefix) );
}
