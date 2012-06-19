#ifndef _FILE_NAME_UTILS_HH_
#define _FILE_NAME_UTILS_HH_
#include <string>

class FileNameUtils {
public:
  static string getSeqFileName( const string& prefix )                { return prefix + "/Sequences"; }
  static string getRoadmapFileName( const string& prefix )            { return prefix + "/Roadmaps"; }
  static string getPregraphFileName( const string& prefix )           { return prefix + "/PreGraph"; }
  static string getInitialGraphFileName( const string& prefix )       { return prefix + "/Graph2"; }
  static string getVelvetLastGraphFileName( const string& prefix )    { return prefix + "/LastGraph"; }
  static string getLastGraphFileName( const string& prefix )          { return prefix + "/meta-velvetg.LastGraph"; }
  static string getScaffoldFileName( const string& prefix )           { return prefix + "/meta-velvetg.contigs.fa"; }
  static string getHistoStatsFileName( const string& prefix )         { return prefix + "/meta-velvetg.Graph2-stats.txt"; } 
  static string getStatsFileName( const string& prefix )              { return prefix + "/meta-velvetg.LastGraph-stats.txt"; }
  static string getSplitStatsFileName( const string& prefix )         { return prefix + "/meta-velvetg.split-stats.txt"; }
  static string getSplitDetailFileName( const string& prefix )        { return prefix + "/meta-velvetg.split-detail.fa"; }
  static string getSubgraphPrefix( const string& prefix )             { return prefix + "/meta-velvetg.subgraph"; }
  static string getAmosFileName( const string& prefix )               { return prefix + "/meta-velvetg.asm.afg"; }
  static string getAlignmentFileName( const string& prefix )          { return prefix + "/meta-velvetg.contig-alignments.psa"; }
  static string getLowCoverageContigsFileName( const string& prefix ) { return prefix + "/meta-velvetg.lowCoverageContigs.fa"; }
  static string getHighCoverageContigsFileName( const string& prefix ){ return prefix + "/meta-velvetg.highCoverageContigs.fa"; }
};

#endif // _FILE_NAME_UTILS_HH_
