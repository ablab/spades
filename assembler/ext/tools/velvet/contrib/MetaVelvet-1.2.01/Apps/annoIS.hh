#ifndef _ANNOIS_HH_
#define _ANNOIS_HH_
#include <stdio.h>
#include <string.h>
#include <string>
#include <iostream>
#include "../VelvetAPI/VelvetHeaders.hh"
#include "../VelvetAPI/VelvetUtils.hh"
#include "../Common/FileUtils.hh"
#include "../ISGraph/ISGraph.hh"
#include "FileNameUtils.hh"
#include "AppsDefs.hh"

#define ANNOIS_DEFAULT_RELIABLE_NODE_LENGTH 500
#define ANNOIS_DEFAULT_FLANKING_LENGTH 100

using namespace std;

class AnnoIS {
  string prefix;
  string nodeFileName;
  string arNodeFileName;
  int reliableNodeLength;
  int flankingLength;
  bool flagUseMV;
  void setDefaultParameters();
  bool checkParameters() const;
  bool checkParameters( int argc, char* argv[] ) const;
  bool checkArgc( int argc ) const;
  bool checkPrefix( int argc, char* argv[] ) const;
  bool checkMandatoryFiles( string prefix ) const;
  bool checkCategory( int cat ) const;
  void printUsage() const;
  bool setParameters( int argc, char* argv[] );
  bool setParameter( char* name, char* val );
  bool setReliableNodeLength( char* val );
  bool setFlankingLength( char* val );
  bool setFlagUseMV( char* val );
  string getGraphFileName( const string& prefix ) const;
public:
  AnnoIS( int argc, char* argv[] );
  vector<ISNode*> loadISNodes() const;
  ISGraph* loadGraph() const;
  void annotateIS( ISGraph* graph, const vector<ISNode*>& isNodes ) const;
  void saveARNodes( const vector<ISNode*>& isNodes ) const;
  string getARNodeFileName() const { return arNodeFileName; }
};

#endif // _ANNOIS_HH_
