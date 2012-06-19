#ifndef _IS_NODE_HH_
#define _IS_NODE_HH_
#include <map>
#include <iostream>
#include "../VelvetAPI/VelvetGraph.hh"
#include "../Utils/Utils.hh"
#include "ISPath.hh"

#define IS_NODE_IO_DELIMITER "\t"
#define IS_NODE_IO_NUM_COLUMNS 7
#define AR_NODE_COLUMN_DELIMITER "|"

using namespace std;

class ISNode {
  IDnum nodeID;
  string isName;
  long offset5, expand5, offset3, expand3, length;
  vector<ISPath*> paths5, paths3;
  vector<IDnum> arNodeIDs5, arNodeIDs3;
  void detectARNodeIDs5();
  void detectARNodeIDs3();
public:
  ISNode( IDnum id, string name, long o5, long e5, long o3, long e3, long l ) 
  : nodeID(id), isName(name), offset5(o5), expand5(e5), offset3(o3), expand3(e3), length(l) {}
  IDnum getNodeID() const { return nodeID; }
  string getISName() const { return isName; }
  long getOffset5() const { return offset5; }
  long getExpand5() const { return expand5; }
  long getOffset3() const { return offset3; }
  long getExpand3() const { return expand3; }
  long getLength() const { return length; }
  void addPath5( ISPath* path ){ paths5.push_back( path ); }
  void addPath3( ISPath* path ){ paths3.push_back( path ); }
  void detectARNodeIDs();
  uint getNumARNodes5() const { return arNodeIDs5.size(); }
  uint getNumARNodes3() const { return arNodeIDs3.size(); }
  string getARNodeColumn5() const;
  string getARNodeColumn3() const;
};

class ISNodeIO {
  static vector<ISNode*> load( ifstream& ifs );
  static bool isValidLine( const string& line );
  static ISNode* parseLine( const string& line );
  static void saveARNodes( ofstream& ofs, const vector<ISNode*>& isNodes );
public:
  static vector<ISNode*> load( const string& filename );
  static void saveARNodes( const string& filename, const vector<ISNode*>& isNodes );
};

#endif // _IS_NODE_HH_
