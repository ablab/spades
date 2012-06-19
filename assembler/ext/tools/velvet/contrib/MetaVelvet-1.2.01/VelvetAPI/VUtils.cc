#include "VUtils.hh"

char* VUtils::getCharFileName( const string& name ){
  char* filename;
  filename = mallocOrExit( name.length()+1, char );
  strcpy( filename, name.c_str() );
  return filename;
}

double VUtils::getNodeDensity( Node* node ){
  Coordinate nodeLength = getNodeLength(node);
  Coordinate nodeCoverage = (getVirtualCoverage(node, 0) + getVirtualCoverage(node, 1));
  return nodeCoverage /(double) nodeLength;
}

int VUtils::getNumOutArcs( Node* node ){
  return simpleArcCount(node);
}

int VUtils::getNumInArcs( Node* node ){
  return simpleArcCount( getTwinNode(node) );
}

vector<Node*> VUtils::getOutNodes( Node* node ){
  vector<Node*> outNodes;
  Arc* tmpArc = getArc( node );
  while( tmpArc != NULL ){
    Node* tmpNode = getDestination(tmpArc);
    outNodes.push_back( tmpNode );
    tmpArc = getNextArc( tmpArc );
  }
  return outNodes;
}

vector<Node*> VUtils::getInNodes( Node* node ){
  vector<Node*> inNodes;
  Arc* tmpArc = getArc( getTwinNode(node) );
  while( tmpArc != NULL ){
    Node* tmpNode = getDestination(tmpArc);
    inNodes.push_back( getTwinNode(tmpNode) );
    tmpArc = getNextArc( tmpArc );
  }
  return inNodes;
}

string VUtils::getTitle( Node* node ){
  string delim = VUTILS_FASTA_TITLE_DELIMITER;
  return Utils::ltoa( getNodeID(node) ) + delim
    + Utils::dtoa( getNodeDensity(node) ) + delim
    + Utils::ltoa( getNodeLength(node) );
}

string VUtils::getSequence( Node* node ){
  string seq = "";
  for( Coordinate i=0 ; i<getNodeLength(node) ; ++i ){
    Nucleotide base = getNucleotideInNode( node, i );
    switch(base){
    case ADENINE:
      seq += "A";
      break;
    case CYTOSINE:
      seq += "C";
      break;
    case GUANINE:
      seq += "G";
      break;
    case THYMINE:
      seq += "T";
      break;
    }
  }
  return seq;
}

double VUtils::getTotalNodeDensity( const vector<Node*>& nodes ){
  double totalDensity = 0.0;
  for( size_t i=0 ; i<nodes.size() ; ++i ){
    totalDensity += getNodeDensity( nodes.at(i) );
  }
  return totalDensity;
}

double VUtils::getMinNodeDensity( const vector<Node*>& nodes ){
  double minDensity = getNodeDensity( nodes.at(0) );
  for( size_t i=1 ; i<nodes.size() ; ++i ){
    double curDensity = getNodeDensity( nodes.at(i) );
    if( curDensity < minDensity ){
      minDensity = curDensity;
    }
  }
  return minDensity;
}

double VUtils::getMaxNodeDensity( const vector<Node*>& nodes ){
  double maxDensity = getNodeDensity( nodes.at(0) );
  for( size_t i=1 ; i<nodes.size() ; ++i ){
    double curDensity = getNodeDensity( nodes.at(i) );
    if( curDensity > maxDensity ){
      maxDensity = curDensity;
    }
  }
  return maxDensity;
}

Node* VUtils::getMinDensityNode( const vector<Node*>& nodes ){
  double minDensity = getNodeDensity( nodes.at(0) );
  Node* minNode     = nodes.at(0);
  for( size_t i=1 ; i<nodes.size() ; ++i ){
    double curDensity = getNodeDensity( nodes.at(i) );
    if( curDensity < minDensity ){
      minDensity = curDensity;
      minNode    = nodes.at(i);
    }
  }
  return minNode;
}

Node* VUtils::getMaxDensityNode( const vector<Node*>& nodes ){
  double maxDensity = getNodeDensity( nodes.at(0) );
  Node* maxNode     = nodes.at(0);
  for( size_t i=1 ; i<nodes.size() ; ++i ){
    double curDensity = getNodeDensity( nodes.at(i) );
    if( curDensity > maxDensity ){
      maxDensity = curDensity;
      maxNode    = nodes.at(i);
    }
  }
  return maxNode;
}

