///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Classes to simplify extracting data from FASTA files.

#ifndef FASTAFILESET_H
#define FASTAFILESET_H

#include "system/MemTracker.h"
#include "String.h"
#include "Vec.h"
#include "VecString.h"
#include "FastaFilestream.h"
#include "FastaNameParser.h"
#include "CompressedSequence.h"


/// Put sequences and names from fasta file into vecbasevector and vecString.
/// This is almost twice as fast as FetchReads (at least when loading
/// many short sequences). It also uses less memory, although still a lot
/// more memory than loading from a fastb file.
/// Ns are accepted and are saved as random bases.
void FastFetchReads(vecbasevector & b, vecString * n, const String &file);

///Convenience method, reads both fastb and fasta files.
inline void LoadReads(vecbasevector & reads, const String & fname) {
  if (fname.Contains(".fastb", -1)) reads.ReadAll(fname);
  else FastFetchReads(reads, 0, fname);
}

/// Put qualities and names from qual file into vecqualvector and vecString.
void FastFetchQuals(vecqualvector & q, vecString * n, const String &file);


/// Simple struct to pair names with data.
/// \class NamedData
template<typename dataT>
struct NamedData {
  const String *mp_name;
  const dataT  *mp_data;
  
  NamedData()
    : mp_name( 0 ), mp_data( 0 ) {}
  
  NamedData( const String *p_name, const dataT *p_data )
    : mp_name( p_name ), mp_data( p_data ) {}
  
  bool operator< ( const NamedData<dataT> &other ) const
  {
    return ( *mp_name < *(other.mp_name) );
  }
};

/// Template for class that knows how to extract a particular kind of
/// data (dataT) from a set of given files and how to store that data
/// in a vector (vecT) and iterate through it (GetNext()) or find
/// entries in it (GetByName()).
/// \class FastaFilesetTemplate

template<typename vecT, typename dataT, typename filestreamT>
class FastaFilesetTemplate
{
 public:
  
  FastaFilesetTemplate( const vec<String>& filenames,
                        FastaNameParser *nameParser = 0,
                        ostream &log = cout );

  /// Returns the number of sequences found.
  int Parse();

  /// Get the number of sequences found.
  int GetSize();
  
  /// Get the number of bases found.
  longlong GetNumberOfBases();
  
  /// Retrieve the next parsed sequence, returning false if the last
  /// sequence has already been returned. Sequences are ordered 
  /// alphabetically by name.
  bool GetNext( String &name, dataT &bases );

  /// Retrieve the sequence at index i, returning false if the index is out
  /// of bounds. Sequences are ordered 
  /// in the order they were in in the original files.
  bool GetAtIndex( typename vecT::size_type i, String &name, dataT &bases );

  /// Retrieve the next parsed sequence's name, returning false if the last
  /// sequence has already been returned. Sequences are ordered 
  /// alphabetically by name. 
  bool GetNextName( String &name );

  /// Start over at beginning of sequences.
  void Reset();
 
  /// If the parsed sequences contain exactly one piece of data with
  /// the given name, copy it to the parameter and return true.
  /// Otherwise, return false.  Has no effect on iteration via
  /// GetNext().
  bool GetByName( const String &name, dataT &bases );

 private:
  longlong EstimateVecSize( longlong totalSeqSize );

  vec<String> m_filenames;
  
  vecString m_names;
  vecT      m_data;

  typedef struct NamedData<dataT> NamedDataT;
  vec<NamedDataT>                    m_sortedData;
  typename vec<NamedDataT>::iterator m_currentDataIter;

  bool m_isParsed;

  FastaNameParser *mp_nameParser;

  ostream &m_log;
};

// Typedefs for commonly used FastaFilesetTemplates.
typedef FastaFilesetTemplate<veccompseq,CompressedSequence,FastaSequenceFilestream> 
FastaSequenceFileset;

typedef FastaFilesetTemplate<vecqualvector,qualvector,FastaQualityFilestream> 
FastaQualityFileset;


/// This class knows how to pair up and iterate through data from a
/// FastaSequenceFileset and a FastaQualityFileset.
class FastaPairedFileset {

 public:
  
  FastaPairedFileset( const vec<String>& sequenceFilenames, 
                      const vec<String>& qualityFilenames,
                      FastaNameParser *nameParser = 0,
                      ostream &log = cout );

  virtual ~FastaPairedFileset() {}

  /// Returns the number of sequences with both sequence and quality data.
  int Parse();

  /// Retrieve the next parsed sequence for which there is both base
  /// and qual data, returning true if such a sequence was found.
  /// Sequences are ordered alphabetically by name.
  bool GetNext( String &name, CompressedSequence &bases, qualvector &quals );
 
  /// If the parsed sequences contain exactly one set each of bases and
  /// quals with the given name, copy them to the parameters and return
  /// true.  Otherwise, return false.  Has no effect on iteration via
  /// GetNext().
  bool GetByName( const String &name, CompressedSequence &bases, qualvector &quals );

  void GetUnmatchedSequenceNames( vec<String> &unmatchedNames );
  void GetUnmatchedQualityNames( vec<String> &unmatchedNames );

 private:
  FastaSequenceFileset m_basesFileset;
  FastaQualityFileset  m_qualsFileset;

  vec<String> m_unmatchedSeqs;
  vec<String> m_unmatchedQuals;

  bool   m_isParsed;

  ostream &m_log;
};

#endif
