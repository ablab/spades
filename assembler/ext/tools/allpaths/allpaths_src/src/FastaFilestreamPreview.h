///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#ifndef FASTAFILESTREAMPREVIEW_H
#define FASTAFILESTREAMPREVIEW_H

#include "Vec.h"

// FastaFilestreamPreview understands just enough about the fasta
// format to be able to count the number of sequences in the specified
// filestream and to know where each sequence starts in that filestream.

// The use of this class by FastaFilestream makes its job much easier
// and faster.

class FastaFilestreamPreview {

public:

  FastaFilestreamPreview(istream& fasta_istream);

  const streampos getMaxSequenceSize();  
  const streampos getStartOffset();  

  vec<streampos>& getSequenceSizes();

private:

  streampos max_sequence_size_, start_offset_;
  vec<streampos> sequence_sizes_;

};

#endif
