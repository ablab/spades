#include "ISPath.hh"

bool ISPath::hasLoop() const {
  map<IDnum, bool> nodeID2flag;
  for( uint i=0 ; i<path.size() ; ++i ){
    IDnum nodeID = path.at(i);
    map<IDnum, bool>::const_iterator it = nodeID2flag.find( nodeID );
    if( it != nodeID2flag.end() ){
      return true;
    }
    nodeID2flag.insert( map<IDnum, bool>::value_type(nodeID, true) );
  }
  return false;
}

string ISPath::getLine() const {
  string line = (string)"";
  for( uint i=0 ; i<path.size() ; ++i ){
    if( i>0 ){ 
      line += "|";
    }
    line += Utils::itoa( path.at(i) );
  }
  line += "\t";
  for( uint i=0 ; i<path.size() ; ++i ){
    if( i ){
      line += "|";
    }
    line += Utils::itoa( nodeLengths.at(i) );
  }
  return line;
}

void ISPath::show() const {
  cerr << "[ISPath] " << getLine() << endl;
}

ISPath* ISPath::instantiateRoot( IDnum nodeID, int nodeLength ){
  ISPath* root = new ISPath();
  root->addNode( nodeID, nodeLength );
  return root;
}

ISPath* ISPath::instantiateNextPath( const ISPath* p, IDnum nextNodeID, long nextNodeLength ){
  ISPath* nextPath = new ISPath();
  for( uint i=0 ; i<p->getPathLength() ; ++i ){
    nextPath->addNode( p->getNodeID(i), p->getNodeLength(i) );
  }
  nextPath->addNode( nextNodeID, nextNodeLength );
  return nextPath;
}
