#include "FastaWriter.hh"

void FastaWriter::open( const string& filename ) {
  cout << prefix << " Writing Fasta file: " << endl
       << "\t" << filename << endl;
  Utils::fileopen( ofs, filename );
}

void FastaWriter::send( const FastaSeq* faSeq ) {
  ++numSequences;
  ofs << '>' << faSeq->getTitle() << endl;
  for( size_t i=0 ; i<faSeq->getLength() ; i+=numBasesPerLine ){
    if( i+numBasesPerLine < faSeq->getLength() ){
      ofs << faSeq->getSequence().substr( i, numBasesPerLine );
    } else {
      ofs << faSeq->getSequence().substr( i );
    }
    ofs << endl;
  }
}

void FastaWriter::close() {
  cout << prefix << " Complete to write Fasta file: " << endl
       << "\t" << numSequences << " sequences were output." << endl;
  ofs.close();
}
