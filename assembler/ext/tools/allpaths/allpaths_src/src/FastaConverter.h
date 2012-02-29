///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#ifndef FASTACONVERTER_H
#define FASTACONVERTER_H

#include "Charvector.h"
#include "CompressedSequence.h"
#include "FastaNameParser.h"
#include "Qualvector.h"

#include <functional>

// FastaConverter understands the fasta format and can extract
// sequence names and sequence data from a null-terminated buffer
// containing fasta data.

template<typename sequenceT>
class FastaConverter {

 public:
  FastaConverter( FastaNameParser* name_parser )
    : name_parser_(name_parser) { }

  virtual ~FastaConverter() {};

  void setNameParser( FastaNameParser* name_parser )
    { 
      name_parser_ = name_parser;
      cout << name_parser_ << endl;
    }

  void extractNameFromBuffer(char* buffer, String &name );

  bool extractDatumFromBuffer(char* buffer, sequenceT &sequence );

  void extractAllFromBuffer(char* buffer, String &name, sequenceT &sequence );

 protected:
  FastaNameParser* name_parser_;
};

// FastaNullConverter knows nothing other than that the data is stored
// in fasta format.

typedef FastaConverter< charvector > FastaNullConverter;


// FastaSequenceConverter specifically understands how to extract
// sequence data stored in fasta format.

typedef FastaConverter< CompressedSequence > FastaSequenceConverter;


// FastaQualityConverter specifically understands how to extract
// quality data stored in fasta format.

typedef FastaConverter< qualvector > FastaQualityConverter;

#endif
