#include "FastaSeq.hh"

FastaSeq* FastaSeq::instantiate( const string& title, const string& sequence ) {
  FastaSeq* seq = new FastaSeq();
  seq->setTitle(title);
  seq->setSequence(sequence);
  return seq;
}

FastaIndex::FastaIndex( const vector<FastaSeq*>& seqs ) {
  for( size_t i=0 ; i<seqs.size() ; ++i ){
    index.insert( map<string, FastaSeq*>::value_type( seqs.at(i)->getTitle(), seqs.at(i) ) );
  }
}

FastaSeq* FastaIndex::getSequence( const string& title ) const {
  map<string, FastaSeq*>::const_iterator it = index.find( title );
  return it->second;
}
