#ifndef _IS_GRAPH_HH_
#define _IS_GRAPH_HH_
#include <math.h>
#include <stack>
#include "../VelvetAPI/VelvetGraph.hh"
#include "ISNode.hh"

using namespace std;

class ISGraph : public VelvetGraph {
  void traverse( ISNode* isNode, int reliableNodeLength );
  void traverse5( ISNode* isNode, int reliableNodeLength );
  void traverse3( ISNode* isNode, int reliableNodeLength );
  bool isEndOfSearch( ISPath* p, int reliableNodeLength );
public:
  ISGraph( const string& seqFileName, const string& graphFileName ) : VelvetGraph( seqFileName, graphFileName ){}
  void annotateIS( const vector<ISNode*>& isNodes, int reliableNodeLength, int flankingLength );
};

#endif // _IS_GRAPH_HH_
