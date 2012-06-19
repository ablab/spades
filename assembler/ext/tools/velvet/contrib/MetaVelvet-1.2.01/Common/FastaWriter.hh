#ifndef _FASTA_WRITER_HH_
#define _FASTA_WRITER_HH_
#include "../Utils/Utils.hh"
#include "FastaSeq.hh"

#define FASTA_WRITER_NUM_BASES_PER_LINE 60

using namespace std;

class FastaWriter {
  string prefix;
  ofstream ofs;
  long numSequences;
  size_t numBasesPerLine;
public:
  FastaWriter( const string& p="" ) : prefix(p), numSequences(0), numBasesPerLine(FASTA_WRITER_NUM_BASES_PER_LINE) {}
  void open( const string& filename );
  void send( const FastaSeq* faSeq );
  void close();
};

#endif // _FASTA_WRITER_HH_
