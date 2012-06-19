#ifndef _V_UTILS_HH_
#define _V_UTILS_HH_
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>
#include "../Utils/Utils.hh"
#include "VelvetHeaders.hh"

#define VUTILS_FASTA_TITLE_DELIMITER "___"

using namespace std;

class VUtils {
public:
  static char* getCharFileName( const string& name );
  static double getNodeDensity( Node* node );
  static int getNumOutArcs( Node* node );
  static int getNumInArcs( Node* node );
  static vector<Node*> getOutNodes( Node* node );
  static vector<Node*> getInNodes( Node* node );
  static string getTitle( Node* node );
  static string getSequence( Node* node );
  static double getTotalNodeDensity( const vector<Node*>& nodes );
  static double getMaxNodeDensity( const vector<Node*>& nodes );
  static double getMinNodeDensity( const vector<Node*>& nodes );
  static Node* getMaxDensityNode( const vector<Node*>& nodes );
  static Node* getMinDensityNode( const vector<Node*>& nodes );
};

#endif // _V_UTILS_HH_
