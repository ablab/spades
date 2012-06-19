#ifndef _META_VELVET_H_HH_
#define _META_VELVET_H_HH_
#include <string>
#include <iostream>
#include "../VelvetAPI/VelvetHeaders.hh"
#include "../VelvetAPI/VelvetUtils.hh"
#include "FileNameUtils.hh"
#include "AppsDefs.hh"

#define META_VELVET_H_FLAG_DOUBLE_STRAND true
#define META_VELVET_H_FLAG_NO_HASH false

using namespace std;

class MetaVelvetH {
  bool flagDoubleStrand;
  bool flagNoHash;
  string prefix;
  int hashLength;
  void setDefaultParameters();
  bool checkParameters( int argc, char* argv[] ) const;
  bool checkArgc( int argc ) const;
  bool checkHashLength( int argc, char* argv[] ) const;
  bool checkMultipleHashLength( int argc, char* argv[] ) const;
  bool checkMaxHashLength( int argc, char* argv[] ) const;
  bool checkMinHashLength( int argc, char* argv[] ) const;
  bool checkEvenHashLength( int argc, char* argv[] ) const;
  void printUsage() const;
  void setParameters( int argc, char* argv[] );
public:
  MetaVelvetH( int argc, char* argv[] );
  void importSequences( int argc, char* argv[] );
  void setRoadmaps();
};

#endif // _META_VELVET_H_HH_
