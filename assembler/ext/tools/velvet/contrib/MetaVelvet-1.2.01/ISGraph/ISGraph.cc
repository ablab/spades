#include "ISGraph.hh"

void ISGraph::annotateIS( const vector<ISNode*>& isNodes, int reliableNodeLength, int flankingLength ){
  for( uint i=0 ; i<isNodes.size() ; ++i ){
    cout << "[ISGraph] " << "Node ID = " << isNodes.at(i)->getNodeID() << endl;
    traverse( isNodes.at(i), reliableNodeLength );
    isNodes.at(i)->detectARNodeIDs();
    cout << "[ISGraph] " << "#. AR Nodes (5' direction) = " << isNodes.at(i)->getNumARNodes5() << endl
	 << "[ISGraph] " << "#. AR Nodes (3' direction) = " << isNodes.at(i)->getNumARNodes3() << endl << endl;
  }
}

void ISGraph::traverse( ISNode* isNode, int reliableNodeLength ){
  traverse5( isNode, reliableNodeLength );
  traverse3( isNode, reliableNodeLength );
}

void ISGraph::traverse5( ISNode* isNode, int reliableNodeLength ){
  stack<ISPath*> dfsStack; 
  ISPath* root = ISPath::instantiateRoot( isNode->getNodeID(), 0 );
  dfsStack.push( root );
  while( dfsStack.size() ){
    ISPath* p  = dfsStack.top(); dfsStack.pop();
    if( isEndOfSearch( p, reliableNodeLength ) ){
      cout << "[ISGraph] " << "Found path (5' direction): " << p->getLine() << endl;
      isNode->addPath5( p );
    } else {
      Node* node = this->getNode( p->getLastNodeID() );
      Node* twin = getTwinNode( node );
      Arc*  arc  = getArc( twin );
      while( arc != NULL ){
	Node* nextNode = getDestination( arc );
	Node* nextTwin = getTwinNode( nextNode );
	ISPath* nextPath = ISPath::instantiateNextPath( p, getNodeID(nextTwin), getNodeLength(nextTwin) );
	if( !nextPath->hasLoop() ){
	  dfsStack.push( nextPath );
	}
	arc = getNextArc( arc );
      }
    }
  }
}

void ISGraph::traverse3( ISNode* isNode, int reliableNodeLength ){
  stack<ISPath*> dfsStack;
  ISPath* root = ISPath::instantiateRoot( isNode->getNodeID(), 0 );
  dfsStack.push( root );
  while( dfsStack.size() ){
    ISPath* p  = dfsStack.top(); dfsStack.pop();
    if( isEndOfSearch(p, reliableNodeLength) ){
      cout << "[ISGraph] " << "Found path (3' direction): " << p->getLine() << endl;
      isNode->addPath3( p );
    } else {
      Node* node = this->getNode( p->getLastNodeID() );
      Arc*  arc  = getArc( node );
      while( arc != NULL ){
	Node* nextNode = getDestination( arc );
	ISPath* nextPath = ISPath::instantiateNextPath( p, getNodeID(nextNode), getNodeLength(nextNode) );
	if( !nextPath->hasLoop() ){
	  dfsStack.push( nextPath );
	}
	arc = getNextArc( arc );
      }
    }
  }
}

bool ISGraph::isEndOfSearch( ISPath* p, int reliableNodeLength ){
  return (p->getLastNodeLength() >= reliableNodeLength);
}

