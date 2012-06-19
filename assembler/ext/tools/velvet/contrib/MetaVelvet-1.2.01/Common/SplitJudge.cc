#include "SplitJudge.hh"

void SplitJudge::setStatsWriter( const string& statsFileName, const string& detailFileName, bool flagReportDetail ){
  this->flagReportDetail = flagReportDetail;
  if( flagReportDetail ){
    SplitDetailWriter* writer = new SplitDetailWriter();
    writer->open( statsFileName, detailFileName );
    statsWriter = writer;
  } else {
    SplitStatsWriter* writer = new SplitStatsWriter();
    writer->open( statsFileName );
    statsWriter = writer;
  }
}

void SplitJudge::setExpectedCoverages( int num, const double* covs ){
  numCoveragePeaks = num;
  for( int i=0 ; i<META_GRAPH_MAX_NUM_COVERAGE_PEAKS ; ++i ){ 
    expectedCoverages[i] = covs[i]; 
  }
}

void SplitJudge::showExpectedCoverages() const {
  for( int i=0 ; i<numCoveragePeaks ; ++i ){
    cout << "[MetaGraph] " << (i+1) << "-th coverage peak = " << expectedCoverages[i] << endl;
  }
}

bool SplitJudge::isSplit( Node* node ) const {
  bool flagSplit = false;
  if( isLocalRepeatStructure( node ) ) {
    SplitStats* stats = new SplitStats( statsWriter->getNextID(), node, flagReportDetail ); 
    if( isPassLengthFilter( node ) ){
      stats->setFlagLength( true );
      if( isRepeatCoverageCondition( node ) ){
	stats->setFlagCov( true );
	if( !flagUseConnections || isConsistentConnections( node, stats ) ){
	  stats->setFlagPE( true );
	  flagSplit = true;
	}
      }
    }
    statsWriter->send( stats );
    delete stats;
  }
  return flagSplit;
}

bool SplitJudge::isLocalRepeatStructure( Node * node ) const {
  if( !(VUtils::getNumOutArcs(node) == 2) ){
    return false;
  }
  if( !(VUtils::getNumInArcs(node) == 2) ){
    return false;
  }
  return true;
}

bool SplitJudge::isPassLengthFilter( Node* node ) const {
  return isPassLengthFilter( VUtils::getInNodes( node ) ) and isPassLengthFilter( VUtils::getOutNodes( node ) );
}

bool SplitJudge::isPassLengthFilter( const vector<Node*>& nodes ) const {
  for( size_t i=0 ; i<nodes.size() ; ++i ){
    Coordinate nodeLen = getNodeLength( nodes.at(i) );
    if( nodeLen < (Coordinate)minSplitLength ){
      return false;
    }
  }
  return true;
}

bool SplitJudge::isRepeatCoverageCondition( Node* node ) const {
  vector<Node*> inNodes  = VUtils::getInNodes( node );
  vector<Node*> outNodes = VUtils::getOutNodes( node );
  double inPeak1  = getNearestPeak( VUtils::getMaxNodeDensity( inNodes ) );
  double inPeak2  = getNearestPeak( VUtils::getMinNodeDensity( inNodes ) );
  double outPeak1 = getNearestPeak( VUtils::getMaxNodeDensity( outNodes ) );
  double outPeak2 = getNearestPeak( VUtils::getMinNodeDensity( outNodes ) );
  if( inPeak1==outPeak1 and inPeak2==outPeak2 and inPeak1>inPeak2 ){
    double aveDensity = ( VUtils::getTotalNodeDensity(inNodes) + VUtils::getTotalNodeDensity(outNodes) ) / 2.0;
    if( VUtils::getNodeDensity(node)<=aveDensity*(1.0+repeatCoverageSD) and VUtils::getNodeDensity(node)>=aveDensity*(1.0-repeatCoverageSD) ){
      return true;
    }
  }
  return false;
}

bool SplitJudge::isConsistentConnections( Node* node, SplitStats* stats ) const {
  size_t numConsistentConnections   = calcNumConsistentConnections( node );
  size_t numInConsistentConnections = calcNumInConsistentConnections( node );
  stats->setNumConsistentConnections( numConsistentConnections );
  stats->setNumInConsistentConnections( numInConsistentConnections );
  return ( numInConsistentConnections <= numNoiseConnections ) and ( numConsistentConnections >= numValidConnections );
}

size_t SplitJudge::calcNumConsistentConnections( Node* node ) const {
  return calcNumConsistentConnectionsForPrimaryNodes(node) + calcNumConsistentConnectionsForNonPrimaryNodes(node);
}

size_t SplitJudge::calcNumConsistentConnectionsForPrimaryNodes( Node* node ) const {
  vector<Node*> inNodes  = VUtils::getInNodes( node );
  vector<Node*> outNodes = VUtils::getOutNodes( node );
  Node* inNode  = VUtils::getMaxDensityNode( inNodes );
  Node* outNode = VUtils::getMaxDensityNode( outNodes );
  return calcNumConnections( inNode, outNode ) + calcNumConnections( getTwinNode(inNode), outNode );
}

size_t SplitJudge::calcNumConsistentConnectionsForNonPrimaryNodes( Node* node ) const {
  vector<Node*> inNodes  = VUtils::getInNodes( node );
  vector<Node*> outNodes = VUtils::getOutNodes( node );
  Node* inNode  = VUtils::getMinDensityNode( inNodes );
  Node* outNode = VUtils::getMinDensityNode( outNodes );
  return calcNumConnections( inNode, outNode ) + calcNumConnections( getTwinNode(inNode), outNode );
}

size_t SplitJudge::calcNumInConsistentConnections( Node* node ) const {
  vector<Node*> inNodes  = VUtils::getInNodes( node );
  vector<Node*> outNodes = VUtils::getOutNodes( node );
  Node* inNode1  = VUtils::getMaxDensityNode( inNodes );
  Node* inNode2  = VUtils::getMinDensityNode( inNodes );
  Node* outNode1 = VUtils::getMaxDensityNode( outNodes );
  Node* outNode2 = VUtils::getMinDensityNode( outNodes );
  return calcNumConnections( inNode1, outNode2 ) + calcNumConnections( getTwinNode(inNode1), outNode2 )
    + calcNumConnections( inNode2, outNode1 ) + calcNumConnections( getTwinNode(inNode2), outNode1 );
}

size_t SplitJudge::calcNumConnections( Node* node1, Node* node2 ) const {
  size_t numLinks = 0;
  for( Connection* link = getConnection(node1) ; link != NULL ; link = getNextConnection(link) ){
    Node* destination = getConnectionDestination( link );
    if( destination == node2 || destination == getTwinNode(node2) ){
      ++numLinks;
    }
  }
  return numLinks;
}

double SplitJudge::getNearestPeak( double coverage ) const {
  double nearestPeak = expectedCoverages[0];
  double nearestDiff = fabs( coverage - nearestPeak );
  for( int i=1 ; i<numCoveragePeaks ; ++i ){
    double curPeak = expectedCoverages[i];
    double curDiff = fabs( coverage - curPeak );
    if( curDiff < nearestDiff ){
      nearestPeak = curPeak;
      nearestDiff = curDiff;
    }
  }
  return nearestPeak;
}

