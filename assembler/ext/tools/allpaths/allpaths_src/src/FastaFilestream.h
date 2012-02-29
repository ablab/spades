///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#ifndef FASTAFILESTREAM_H
#define FASTAFILESTREAM_H

#include "FastaFilestreamPreview.h"
#include "FastaConverter.h"
#include "FastaVerifier.h"
#include "Charvector.h"
#include "VecString.h"

#include <functional>

/// FastaFilestream can parse the fasta file specified in the ctor.

// superclass

template<typename vecT, typename seqT, typename converterT, typename verifierT>
class FastaFilestream {

 public:
  FastaFilestream( const String& filename,
		   FastaNameParser* name_parser );

  FastaFilestream( const FastaFilestream& original );
  FastaFilestream& operator= ( const FastaFilestream& original );

  virtual ~FastaFilestream() 
  {
    if (preview_ptr_) 
      delete preview_ptr_; 
  }

  /// Get estimated number of sequences.
  int estimatedSize();

  /// Get estimated amount of data in bases.
  longlong estimatedData();


  /// Verify the format of the file specified in the constructor.
  bool verify();

  /// Get the names of the sequences, but Reserve first!
  /// This routine does not try to resize names, but just uses
  /// push_back, so if you don't Reserve you will get quadratic performance.
  void getOnlyNames( vecString &names ) 
  { this->parse_( names, 0, 0 ); }
    
  /// Gets all the sequences in the file, but Reserve first!.
  /// This routine does not try to resize names or sequences, but just uses
  /// push_back, so if you don't Reserve you will get quadratic performance.
  void parse( vecString &names, vecT &sequences ) 
  { this->parse_( names, &sequences, 0 ); }

  /// Gets the sequences with the specified indices, but Reserve first!
  /// Reads sequences from the file specified in the constructor.  
  /// The indices vector must be sorted.
  /// This routine does not try to resize names or sequences, but just uses
  /// push_back, so if you don't Reserve you will get quadratic performance.
  void parseSubset( const vec<int>& indices, vecString &names, vecT &sequences )
  { this->parse_( names, &sequences, &indices ); }

 private:
  String filename_;
  FastaFilestreamPreview* preview_ptr_;
  bool parsed_;
  bool needs_pipe_;

  FastaNameParser* name_parser_;
  
  void preview_( istream* fasta_istream );
  
  /// This method is very slow (quadratic) unless  names and p_sequences
  /// have been appropriately Reserve()d first.
  /// If p_sequences is 0, don't save the sequences.
  /// If p_indices is 0, get all data, else just grab the sequences specified.
  void parse_( vecString &names, vecT *p_sequences, const vec<int> *p_indices );

  istream* getIstream_() const;
  void resetIstream_( istream*& fasta_istream ) const;
  void closeIstream_( istream* fasta_istream ) const;

  const String getBasename_( const String& filename ) const;
};


typedef 
FastaFilestream<veccharvector,charvector,FastaNullConverter,FastaNullVerifier>
FastaRawFilestream;

typedef
FastaFilestream<veccompseq,CompressedSequence,FastaSequenceConverter,FastaSequenceVerifier>
FastaSequenceFilestream;

typedef
FastaFilestream<vecqualvector,qualvector,FastaQualityConverter,FastaQualityVerifier>
FastaQualityFilestream;



template <typename filestreamT>
class FastaFilestreamBuilder : unary_function<String, filestreamT>
{
 public:
  FastaFilestreamBuilder( FastaNameParser* name_parser )
    : name_parser_(name_parser) { }
  
  filestreamT operator() (const String& filename)
  { 
    return filestreamT( filename, name_parser_ ); 
  }
  
 private:
  FastaNameParser* name_parser_;
};


typedef FastaFilestreamBuilder<FastaRawFilestream>      FastaRawFilestreamBuilder;
typedef FastaFilestreamBuilder<FastaSequenceFilestream> FastaSequenceFilestreamBuilder;
typedef FastaFilestreamBuilder<FastaQualityFilestream>  FastaQualityFilestreamBuilder;


#endif
