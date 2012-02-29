///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "FastaFileset.h"

#include <strstream>
#include <algorithm>
#include <numeric>
#include <utility>
#include <sys/types.h>
#include <dirent.h>

template<typename vecT, typename dataT, typename streamT>
FastaFilesetTemplate<vecT,dataT,streamT>::FastaFilesetTemplate( const vec<String>& filenames,
                                                                FastaNameParser *nameParser,
                                                                ostream &log )
  : m_filenames( filenames ),
    m_isParsed( false ),
    mp_nameParser( nameParser ),
    m_log( log )
{
  if ( mp_nameParser == 0 )
      mp_nameParser = new LastWordParser;
}

template<typename vecT, typename dataT, typename streamT>
int
FastaFilesetTemplate<vecT,dataT,streamT>::GetSize( )
{
  if ( ! m_isParsed ) this->Parse();

  return m_sortedData.size();
}

template<typename vecT, typename dataT, typename streamT>
longlong
FastaFilesetTemplate<vecT,dataT,streamT>::GetNumberOfBases( )
{
  if ( ! m_isParsed ) this->Parse();

  longlong nBases = 0;
  for ( unsigned int dataIdx = 0; dataIdx < m_sortedData.size(); ++dataIdx )
    nBases += m_sortedData[ dataIdx ].mp_data->size();

  return nBases;
}

template<typename vecT, typename dataT, typename streamT>
bool
FastaFilesetTemplate<vecT,dataT,streamT>::GetNext( String &name, dataT &datum )
{
  if ( ! m_isParsed ) this->Parse();

  if ( m_currentDataIter == m_sortedData.end() )
    return false;

  name  = *(m_currentDataIter->mp_name);
  datum = *(m_currentDataIter->mp_data);
  ++m_currentDataIter;

  return true;
}

// TODO: Potentially dangerous truncation of ID
template<typename vecT, typename dataT, typename streamT>
bool
FastaFilesetTemplate<vecT,dataT,streamT>::GetAtIndex( typename vecT::size_type i,
						      String &name,
						      dataT &datum )
{
  if ( ! m_isParsed ) this->Parse();
  if ( i >= m_data.size() ) {
    return false;
  }

  name  = m_names[i];
  datum = m_data[i];
  return true;
}

template<typename vecT, typename dataT, typename streamT>
bool
FastaFilesetTemplate<vecT,dataT,streamT>::GetNextName( String &name )
{
  if ( ! m_isParsed ) this->Parse();

  if ( m_currentDataIter == m_sortedData.end() )
    return false;

  name  = *(m_currentDataIter->mp_name);
  ++m_currentDataIter;

  return true;
}

template<typename vecT, typename dataT, typename streamT>
void
FastaFilesetTemplate<vecT,dataT,streamT>::Reset()
{
  if ( ! m_isParsed ) this->Parse();

  m_currentDataIter = m_sortedData.begin();
}

template<typename vecT, typename dataT, typename streamT>
bool
FastaFilesetTemplate<vecT,dataT,streamT>::GetByName( const String &name, dataT &datum )
{
  if ( ! m_isParsed ) this->Parse();

  pair<typename vec<NamedDataT>::iterator,typename vec<NamedDataT>::iterator> dataRange;
  dataRange = equal_range( m_sortedData.begin(), m_sortedData.end(),
                           NamedData<dataT>(&name,0) );

  if ( distance( dataRange.first, dataRange.second ) != 1 )
    return false;

  datum = *(dataRange.first->mp_data);
  return true;
}

template<>
longlong
FastaFilesetTemplate<veccompseq,CompressedSequence,FastaSequenceFilestream>::EstimateVecSize( longlong totalSeqSize )
{
  // average of 1 char per value, 5 values per unit of storage in veccompseq
  return totalSeqSize/5;
}

template<>
longlong
FastaFilesetTemplate<vecqualvector,qualvector,FastaQualityFilestream>::EstimateVecSize( longlong totalSeqSize )
{
  // average of 2.5 chars per value, 1 value per unit of storage in vecqualvector
  return totalSeqSize*2/5;
}


template<typename vecT, typename dataT, typename streamT>
int
FastaFilesetTemplate<vecT,dataT,streamT>::Parse()
{
  m_log << Date() << ": Scanning FASTA files." << endl;

  // read in sequences
  vec<streamT> filestreams;
  transform( m_filenames.begin(), m_filenames.end(),
	     back_inserter(filestreams),
	     FastaFilestreamBuilder<streamT>(mp_nameParser) );

  unsigned int totalSequences = 0;
  for ( unsigned int ii = 0; ii < filestreams.size(); ++ii )
    totalSequences += filestreams[ii].estimatedSize();

  longlong totalSeqData = 0;
  for ( unsigned int ii = 0; ii < filestreams.size(); ++ii )
    totalSeqData += filestreams[ii].estimatedData();

  m_log << Date() << ": Preparing to read " << totalSequences << " FASTA entries." << endl;

  // arbitrary overly-large constant :)
  const longlong averageNameSize = 50;
  m_data.Reserve( this->EstimateVecSize( totalSeqData ), totalSequences );
  m_names.Reserve( totalSequences * averageNameSize, totalSequences );

  for ( unsigned int ii = 0; ii < filestreams.size(); ++ii )
    filestreams[ii].parse( m_names, m_data );

  // Set up vectors of objects that tie data to names.
  m_sortedData.resize( m_names.size() );
  for ( size_t seqIdx = 0; seqIdx < m_names.size(); ++seqIdx )
  {
    m_sortedData[ seqIdx ].mp_name = &m_names[seqIdx];
    m_sortedData[ seqIdx ].mp_data = &m_data[seqIdx];
  }
  sort( m_sortedData.begin(), m_sortedData.end() );

  m_currentDataIter = m_sortedData.begin();

  m_log << Date() << ": Done." << endl;

  m_isParsed = true;

  PRINT_TO( m_log, m_data.sumSizes() );

  return m_sortedData.size();
}

#define INSTANTIATE(vecT,dataT,streamT) \
template \
FastaFilesetTemplate<vecT,dataT,streamT>::FastaFilesetTemplate( const vec<String>& filenames, \
                                                                FastaNameParser *nameParser, \
                                                                ostream &log ); \
template int FastaFilesetTemplate<vecT,dataT,streamT>::GetSize( ); \
template longlong FastaFilesetTemplate<vecT,dataT,streamT>::GetNumberOfBases( ); \
template bool FastaFilesetTemplate<vecT,dataT,streamT>::GetNext( String&, dataT& ); \
template bool FastaFilesetTemplate<vecT,dataT,streamT>::GetAtIndex( vecT::size_type, String&, dataT& ); \
template bool FastaFilesetTemplate<vecT,dataT,streamT>::GetNextName( String& ); \
template void FastaFilesetTemplate<vecT,dataT,streamT>::Reset(); \
template bool FastaFilesetTemplate<vecT,dataT,streamT>::GetByName( const String&, dataT& ); \
template int FastaFilesetTemplate<vecT,dataT,streamT>::Parse();

INSTANTIATE( veccompseq,CompressedSequence,FastaSequenceFilestream )
INSTANTIATE( vecqualvector,qualvector,FastaQualityFilestream )



FastaPairedFileset::FastaPairedFileset( const vec<String>& sequenceFilenames,
                                        const vec<String>& qualityFilenames,
                                        FastaNameParser *nameParser,
                                        ostream &log )
  : m_basesFileset( sequenceFilenames, nameParser, log ),
    m_qualsFileset( qualityFilenames, nameParser, log ),
    m_isParsed( false ),
    m_log( log )
{
}


bool
FastaPairedFileset::GetNext( String &name, CompressedSequence &bases, qualvector &quals )
{
  // Tie basesName to name, so we don't need to copy it if we're successful.
  String &basesName = name;
  String qualsName;

  if ( ! m_isParsed ) this->Parse();

  bool basesFound = m_basesFileset.GetNext( basesName, bases );
  bool qualsFound = m_qualsFileset.GetNext( qualsName, quals );

  while ( basesFound && qualsFound )
  {
    if ( basesName < qualsName )
      basesFound = m_basesFileset.GetNext( basesName, bases );
    else if ( qualsName < basesName )
      qualsFound = m_qualsFileset.GetNext( qualsName, quals );
    else
      break;
  }

  return ( basesFound && qualsFound );
}


bool
FastaPairedFileset::GetByName( const String &name, CompressedSequence &bases, qualvector &quals )
{
  if ( ! m_isParsed ) this->Parse();

  return ( m_basesFileset.GetByName( name, bases ) &&
           m_qualsFileset.GetByName( name, quals ) );
}


void
FastaPairedFileset::GetUnmatchedSequenceNames( vec<String> &unmatchedNames )
{
  if ( ! m_isParsed ) this->Parse();
  unmatchedNames = m_unmatchedSeqs;
}

void
FastaPairedFileset::GetUnmatchedQualityNames( vec<String> &unmatchedNames )
{
  if ( ! m_isParsed ) this->Parse();
  unmatchedNames = m_unmatchedQuals;
}


int
FastaPairedFileset::Parse()
{
  // Read data.
  int numSeqsParsed = m_basesFileset.Parse();
  int numQualsParsed = m_qualsFileset.Parse();

  // Find names with missing data.
  m_log << Date() << ": Pairing sequences with quality scores." << endl;

  String basesName, qualsName;

  bool basesFound = m_basesFileset.GetNextName( basesName );
  bool qualsFound = m_qualsFileset.GetNextName( qualsName );

  // Walk through both filesets, saving the names of sequences that
  // don't have quality scores and of quality scores that don't have
  // sequences.
  while ( basesFound && qualsFound )
  {
    if ( basesName < qualsName )
    {
      m_unmatchedSeqs.push_back( basesName );
      basesFound = m_basesFileset.GetNextName( basesName );
    }
    else if ( qualsName < basesName )
    {
      m_unmatchedQuals.push_back( qualsName );
      qualsFound = m_qualsFileset.GetNextName( qualsName );
    }
    else
    {
      basesFound = m_basesFileset.GetNextName( basesName );
      qualsFound = m_qualsFileset.GetNextName( qualsName );
    }
  }

  while ( basesFound )
  {
    m_unmatchedSeqs.push_back( basesName );
    basesFound = m_basesFileset.GetNextName( basesName );
  }

  while ( qualsFound )
  {
    m_unmatchedQuals.push_back( qualsName );
    qualsFound = m_qualsFileset.GetNextName( qualsName );
  }

  // Sanity check.
  ForceAssertEq( numSeqsParsed - m_unmatchedSeqs.size(),
                 numQualsParsed - m_unmatchedQuals.size() );

  // Set iterators for use by GetNext().
  m_basesFileset.Reset();
  m_qualsFileset.Reset();

  m_log << Date() << ": Done." << endl;

  m_isParsed = true;

  // Return number of sequences with complete data.
  return numSeqsParsed - m_unmatchedSeqs.size();
}

void FastFetchReads(vecbasevector & b, vecString * n, const String &file) {
  vec<String> fnames;
  fnames.push_back(file);
  FullNameParser name_parser;
  Ofstream(nolog, "/dev/null");
  FastaSequenceFileset filestream( fnames, &name_parser, nolog );
  filestream.Parse();
  b.clear();
  b.Reserve(filestream.GetNumberOfBases()/4 + 100, filestream.GetSize());
  if (n) n->clear();

  // now put the data into b and n.
  CompressedSequence compseq;
  String name;
  basevector tempb;
  filestream.Reset();
  const int S = filestream.GetSize();
  for (int i=0; i != S; ++i) {
    filestream.GetAtIndex(i,name,compseq);
    compseq.asBasevector(tempb);
    b.push_back(tempb);
    if (n) n->push_back_reserve(name);
  }
}


void FastFetchQuals(vecqualvector & q, vecString * n, const String &file) {
  vec<String> fnames;
  fnames.push_back(file);
  FullNameParser name_parser;
  Ofstream(nolog, "/dev/null");
  FastaQualityFileset filestream( fnames, &name_parser, nolog );
  filestream.Parse();
  q.clear();
  q.Reserve(filestream.GetNumberOfBases(), filestream.GetSize());
  if (n) n->clear();

  // now put the data into q and n.
  String name;
  qualvector tempq;
  filestream.Reset();
  const int S = filestream.GetSize();
  for (int i=0; i != S; ++i) {
    filestream.GetAtIndex(i,name,tempq);
    q.push_back(tempq);
    if (n) n->push_back_reserve(name);
  }
}
