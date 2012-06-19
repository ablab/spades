#ifndef _FASTA_SEQ_HH_
#define _FASTA_SEQ_HH_
#include <string>
#include <vector>
#include <map>

using namespace std;

class FastaSeq {
  string title;
  string sequence;
public:
  void setTitle( const string& t ){ title = t; }
  void setSequence( const string& s ){ sequence = s; }
  void appendSequence( const string& s ){ sequence += s; }
  string getTitle() const { return title; }
  string getSequence() const { return sequence; }
  size_t getLength() const { return sequence.length(); }
  static FastaSeq* instantiate( const string& title, const string& sequence );
};

class FastaIndex {
  map<string, FastaSeq*> index;
public:
  FastaIndex( const vector<FastaSeq*>& seqs );
  FastaSeq* getSequence( const string& title ) const;
};

#endif // _FASTA_SEQ_HH_
