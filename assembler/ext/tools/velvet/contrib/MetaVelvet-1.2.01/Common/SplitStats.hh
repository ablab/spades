#ifndef _SPLIT_STATS_HH_
#define _SPLIT_STATS_HH_
#include "../Utils/Utils.hh"
#include "../VelvetAPI/VUtils.hh"
#include "../VelvetAPI/VelvetGraph.hh"
#include "FastaWriter.hh"

#define SPLIT_STATS_DELIMITER "\t"
#define SPLIT_STATS_TITLE_DELIMITER "___"

using namespace std;

class SplitStats {
protected:
  size_t splitID;
  IDnum repID, pinID, poutID, sinID, soutID;
  double repCov, pinCov, poutCov, sinCov, soutCov;
  Coordinate repLen, pinLen, poutLen, sinLen, soutLen;
  size_t numConsistentConnections, numInConsistentConnections;
  bool flagLen, flagCov, flagPE;
  FastaSeq *repSeq, *pinSeq, *poutSeq, *sinSeq, *soutSeq;
  FastaSeq* node2fasta( Node* node, const string& status );
  string node2title( Node* node, const string& status );
public:
  SplitStats( size_t splitID, Node* node, bool flagReportDetail );
  void setFlagLength( bool flag ){ flagLen = flag; }
  void setFlagCov( bool flag ) { flagCov = flag; }
  void setFlagPE( bool flag ) { flagPE = flag; }
  void setNumConsistentConnections( size_t num ) { numConsistentConnections = num; }
  void setNumInConsistentConnections( size_t num ) { numInConsistentConnections = num; }
  string getLine() const;
  static string getHeaderLine();
  void sendSequences( FastaWriter* writer ) const;
  ~SplitStats(){
    if( repSeq )  delete repSeq;
    if( pinSeq )  delete pinSeq;
    if( poutSeq ) delete poutSeq;
    if( sinSeq )  delete sinSeq;
    if( soutSeq ) delete soutSeq;
  }
};


class SplitStatsWriter {
protected:
  size_t numStats;
  ofstream ofs;
public:
  void open( const string& filename );
  virtual void send( const SplitStats* stats );
  virtual void close() { ofs.close(); }
  size_t getNextID() const;
};

class SplitDetailWriter : public SplitStatsWriter {
  FastaWriter* detailWriter;
public:
  void open( const string& statsFileName, const string& detailFileName );
  void send( const SplitStats* stats );
  void close() { ofs.close(); detailWriter->close(); delete detailWriter; }
};

#endif // _SPLIT_STATS_HH_
