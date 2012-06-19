#include "SplitStats.hh"

SplitStats::SplitStats( size_t splitID, Node* repNode, bool flagReportDetail ) {
  this->splitID = splitID;
  vector<Node*> inNodes  = VUtils::getInNodes( repNode );
  vector<Node*> outNodes = VUtils::getOutNodes( repNode );
  Node* pinNode  = VUtils::getMaxDensityNode( inNodes );
  Node* poutNode = VUtils::getMaxDensityNode( outNodes );
  Node* sinNode  = VUtils::getMinDensityNode( inNodes );
  Node* soutNode = VUtils::getMinDensityNode( outNodes );
  repID   = getNodeID( repNode ); 
  pinID   = getNodeID( pinNode );
  poutID  = getNodeID( poutNode );
  sinID   = getNodeID( sinNode );
  soutID  = getNodeID( soutNode );
  repCov  = VUtils::getNodeDensity( repNode );
  pinCov  = VUtils::getNodeDensity( pinNode );
  poutCov = VUtils::getNodeDensity( poutNode );
  sinCov  = VUtils::getNodeDensity( sinNode );
  soutCov = VUtils::getNodeDensity( soutNode );
  repLen  = getNodeLength( repNode );
  pinLen  = getNodeLength( pinNode );
  poutLen = getNodeLength( poutNode );
  sinLen  = getNodeLength( sinNode );
  soutLen = getNodeLength( soutNode );
  numConsistentConnections   = 0;
  numInConsistentConnections = 0;
  flagLen = false;
  flagCov = false;
  flagPE  = false;
  if( flagReportDetail ){
    repSeq  = node2fasta( repNode,  "RepeatNode" );
    pinSeq  = node2fasta( pinNode,  "PrimaryInNode" );
    poutSeq = node2fasta( poutNode, "PrimaryOutNode" );
    sinSeq  = node2fasta( sinNode,  "SecondaryInNode" );
    soutSeq = node2fasta( soutNode, "SecondaryOutNode" );
  } else {
    repSeq  = NULL;
    pinSeq  = NULL;
    poutSeq = NULL;
    sinSeq  = NULL;
    soutSeq = NULL;
  }
}

FastaSeq* SplitStats::node2fasta( Node* node, const string& status ){
  return FastaSeq::instantiate( node2title(node, status), VUtils::getSequence(node) );
}

string SplitStats::node2title( Node* node, const string& status ){
  string delim   = SPLIT_STATS_TITLE_DELIMITER;
  size_t splitID = this->splitID;
  IDnum  nodeID  = getNodeID( node );
  double cov     = VUtils::getNodeDensity( node );
  Coordinate len = getNodeLength( node );
  return Utils::itoa(splitID) + delim + status + delim + Utils::ltoa(nodeID) + delim + Utils::dtoa(cov) + delim + Utils::ltoa(len);
}

string SplitStats::getHeaderLine() {
  string delim = SPLIT_STATS_DELIMITER;
  return "split_id" + delim
    + "rep_id" + delim
    + "pin_id"   + delim
    + "pout_id"  + delim
    + "sin_id"   + delim
    + "sout_id"  + delim
    + "rep_cov"  + delim
    + "pin_cov"  + delim
    + "pout_cov" + delim
    + "sin_cov"  + delim
    + "sout_cov" + delim
    + "rep_len"  + delim
    + "pin_len"  + delim
    + "pout_len" + delim
    + "sin_len"  + delim
    + "sout_len" + delim
    + "num_consistent_connections"   + delim
    + "num_inconsistent_connections" + delim
    + "flag_pass_length_filter"      + delim
    + "flag_pass_cov_filter"         + delim
    + "flag_pass_pe_filter";
}

string SplitStats::getLine() const {
  string delim = SPLIT_STATS_DELIMITER;
  return Utils::itoa(splitID) + delim
    + Utils::ltoa(repID)    + delim
    + Utils::ltoa(pinID)    + delim
    + Utils::ltoa(poutID)   + delim
    + Utils::ltoa(sinID)    + delim
    + Utils::ltoa(soutID)   + delim
    + Utils::dtoa(repCov)   + delim
    + Utils::dtoa(pinCov)   + delim
    + Utils::dtoa(poutCov)  + delim
    + Utils::dtoa(sinCov)   + delim
    + Utils::dtoa(soutCov)  + delim
    + Utils::ltoa(repLen)   + delim
    + Utils::ltoa(pinLen)   + delim
    + Utils::ltoa(poutLen)  + delim
    + Utils::ltoa(sinLen)   + delim
    + Utils::ltoa(soutLen)  + delim
    + Utils::itoa(numConsistentConnections)   + delim
    + Utils::itoa(numInConsistentConnections) + delim
    + Utils::itoa(flagLen)  + delim
    + Utils::itoa(flagCov)  + delim
    + Utils::itoa(flagPE);
}

void SplitStats::sendSequences( FastaWriter* writer ) const {
  writer->send( repSeq );
  writer->send( pinSeq );
  writer->send( poutSeq );
  writer->send( sinSeq );
  writer->send( soutSeq );
}

void SplitStatsWriter::open( const string& filename ) {
  numStats = 0;
  Utils::fileopen( ofs, filename );
  ofs << SplitStats::getHeaderLine() << endl;
}

void SplitStatsWriter::send( const SplitStats* stats ) {
  ofs << stats->getLine() << endl;
  ++numStats;
}

void SplitDetailWriter::open( const string& statsFileName, const string& detailFileName ){
  numStats = 0;
  Utils::fileopen( ofs, statsFileName );
  ofs << SplitStats::getHeaderLine() << endl;
  detailWriter = new FastaWriter();
  detailWriter->open( detailFileName );
}

void SplitDetailWriter::send( const SplitStats* stats ){
  ofs << stats->getLine() << endl;
  stats->sendSequences( detailWriter );
  ++numStats;
}

size_t SplitStatsWriter::getNextID() const {
  return numStats + 1;
}
