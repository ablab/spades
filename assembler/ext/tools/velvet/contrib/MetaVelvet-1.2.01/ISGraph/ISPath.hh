#ifndef _IS_PATH_HH_
#define _IS_PATH_HH_
#include <map>
#include "../VelvetAPI/VelvetGraph.hh"
#include "../Utils/Utils.hh"

using namespace std;

class ISPath {
  vector<IDnum> path;
  vector<int> nodeLengths;
public:
  void addNode( IDnum nodeID, int nodeLength ){ path.push_back(nodeID); nodeLengths.push_back(nodeLength); }
  IDnum getNodeID( uint i ) const { return path.at(i); }
  int getNodeLength( uint i ) const { return nodeLengths.at(i); }
  IDnum getLastNodeID() const { return path.at( path.size()-1 ); }
  int getLastNodeLength() const { return nodeLengths.at( path.size()-1 ); }
  uint getPathLength() const { return path.size(); }
  bool hasLoop() const;
  string getLine() const;
  void show() const;
  static ISPath* instantiateRoot( IDnum nodeID, int nodeLength );
  static ISPath* instantiateNextPath( const ISPath* prevPath, IDnum nodeID, long nodeLength );
};



#endif // _IS_PATH_HH_
