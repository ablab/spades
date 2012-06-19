#ifndef _VELVET_UTILS_HH_
#define _VELVET_UTILS_HH_
#include "VUtils.hh"

using namespace std;

class VelvetUtils {
public:
  static void importSequences( const string& seqFileName, int argc, char* argv[], boolean& flagDoubleStrand, boolean& flagNoHash );
  static void setRoadmaps( const string& seqFileName, const string& roadmapFileName, int hashLength, boolean flagDoubleStrand, boolean noHash );
  static void constructPregraph( const string& seqFileName, const string& roadmapFileName, const string& pregraphFileName );
};

#endif // _VELVET_UTILS_HH_
