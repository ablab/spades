#include "ISNode.hh"

void ISNode::detectARNodeIDs(){
  detectARNodeIDs5();
  detectARNodeIDs3();
}

void ISNode::detectARNodeIDs5(){
  map<IDnum, bool> nodeID2flag;
  for( uint i=0 ; i<paths5.size() ; ++i ){
    IDnum nodeID = paths5.at(i)->getLastNodeID();
    nodeID2flag.insert( map<IDnum, bool>::value_type(nodeID, true) );
  }
  for( map<IDnum, bool>::const_iterator it=nodeID2flag.begin() ; it!=nodeID2flag.end() ; ++it ){
    IDnum nodeID = it->first;
    arNodeIDs5.push_back( nodeID );
  }
}

void ISNode::detectARNodeIDs3(){
  map<IDnum, bool> nodeID2flag;
  for( uint i=0 ; i<paths3.size() ; ++i ){
    IDnum nodeID = paths3.at(i)->getLastNodeID();
    nodeID2flag.insert( map<IDnum, bool>::value_type(nodeID, true) );
  }
  for( map<IDnum, bool>::const_iterator it=nodeID2flag.begin() ; it!=nodeID2flag.end() ; ++it ){
    IDnum nodeID = it->first;
    arNodeIDs3.push_back( nodeID );
  }
}

string ISNode::getARNodeColumn5() const {
  string col;
  for( uint i=0 ; i<arNodeIDs5.size() ; ++i ){
    if(i){ col += AR_NODE_COLUMN_DELIMITER; }
    col += Utils::itoa( arNodeIDs5.at(i) );
  }
  return col;
}

string ISNode::getARNodeColumn3() const {
  string col;
  for( uint i=0 ; i<arNodeIDs3.size() ; ++i ){
    if(i){ col += AR_NODE_COLUMN_DELIMITER; }
    col += Utils::itoa( arNodeIDs3.at(i) );
  }
  return col;
}

vector<ISNode*> ISNodeIO::load( const string& filename ){
  ifstream ifs;
  Utils::fileopen( ifs, filename, "[ISNodeIO] " );
  return load(ifs);
}

vector<ISNode*> ISNodeIO::load( ifstream& ifs ){
  vector<ISNode*> isNodes;
  string line;
  while( getline(ifs, line) )
    if( isValidLine(line) )
      isNodes.push_back( parseLine(line) );
  ifs.close();
  return isNodes;
}

bool ISNodeIO::isValidLine( const string& line ){
  return ( Utils::split(line, IS_NODE_IO_DELIMITER).size() == IS_NODE_IO_NUM_COLUMNS );
}

ISNode* ISNodeIO::parseLine( const string& line ){
  vector<string> cols = Utils::split( line, IS_NODE_IO_DELIMITER );
  IDnum id     = atoi( cols.at(0).c_str() );
  string name  = cols.at(1);
  long offset5 = atoi( cols.at(2).c_str() );
  long expand5 = atoi( cols.at(3).c_str() );
  long offset3 = atoi( cols.at(4).c_str() );
  long expand3 = atoi( cols.at(5).c_str() );
  long length  = atoi( cols.at(6).c_str() );
  return new ISNode( id, name, offset5, expand5, offset3, expand3, length );
}

void ISNodeIO::saveARNodes( const string& filename, const vector<ISNode*>& isNodes ){
  ofstream ofs;
  Utils::fileopen( ofs, filename, "[ISNodeIO] " );
  return saveARNodes( ofs, isNodes );
}

void ISNodeIO::saveARNodes( ofstream& ofs, const vector<ISNode*>& isNodes ){
  for( uint i=0 ; i<isNodes.size() ; ++i ){
    const ISNode* isNode = isNodes.at(i);
    ofs << isNode->getNodeID() << "\t" 
	<< isNode->getISName() << "\t"
	<< isNode->getARNodeColumn5() << "\t"
	<< isNode->getARNodeColumn3() << endl;
  }
}
