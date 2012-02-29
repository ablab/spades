/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef PATHS_PDF_ENTRY
#define PATHS_PDF_ENTRY

#include "CommonSemanticTypes.h"
#include "String.h"
#include "Vec.h"
#include "feudal/MasterVec.h"
#include "feudal/SerfVec.h"
#include <ostream>

// Struct: pdf_entry
//
// The struct pdf_entry is meant to replace pair<copy_num_t,prob_t> as a way
// of holding the predicted copy number of unipaths (for a given copy
// number the probability that the unipath has that copy number).  
// Its natural habitat is the VecPdfEntryVec.
//
// pair<copy_num_t,prob_t> itself is bad because alignment forces it to have 
// padding on a 64-bit machine but not on a 32-bit, so a double vector of these
// stored on disk is not cross-compatible.
// TODO: consider halving the size of this structure by storing the probability
// as a float (which would also let us get rid of the padding).
struct pdf_entry {
  // field names are meant to minimize code changes from pair<int,double>
public:
  copy_num_t first;   // int
private:
  int PADDING_NEEDED_BY_64_BIT_MACHINES__DO_NOT_USE;
public:
  prob_t second;   // double

  pdf_entry() {}
  pdf_entry(copy_num_t i, prob_t d) : first(i), second(d) { }

  // No one will ever use these,  but just for show:
  copy_num_t NumCopies() const { return first; }
  prob_t Prob() const { return second; }

  friend Bool operator==( const pdf_entry& p1, const pdf_entry& p2 )
  {    return p1.first == p2.first && p1.second == p2.second;    }
};

TRIVIALLY_SERIALIZABLE(pdf_entry);
typedef SerfVec<pdf_entry> PdfEntryVec;
typedef MasterVec<PdfEntryVec> VecPdfEntryVec;

inline std::ostream& operator<<( std::ostream& s, PdfEntryVec const& v )
{
    s << "[";
    PdfEntryVec::const_iterator end(v.end());
    for ( PdfEntryVec::const_iterator itr(v.begin()); itr != end; ++itr )
        s << " " << itr->first << ":" << itr->second;
    s << " ]";
    return s;
}

// Function to find the most likely value in a set of pdf_entries
inline void
GetMostLikelyValue( int & value, const PdfEntryVec& copy_numbers )
{
  value = -1;
  double prob = 0;
  
  // Choose the value with the highest probability
  for ( PdfEntryVec::size_type j = 0; j < copy_numbers.size( ); j++ ) {
    if ( copy_numbers[j].second > prob ) {
      value = copy_numbers[j].first;
      prob = copy_numbers[j].second;
    }
  }
}

// Wrapper around GetMostLikelyValue.
inline int EstimatedCN( const PdfEntryVec& copy_numbers )
{
  int cn = -1;
  GetMostLikelyValue( cn, copy_numbers );
  return cn;
}

// Load CN values and probabilities.
inline void
LoadCopyNumbers( const String &in_file, vec<int> &CNval )
{
  CNval.clear( );

  VecPdfEntryVec pdfs( in_file.c_str() );
  
  CNval.resize( pdfs.size( ), -1 );
  for (VecPdfEntryVec::size_type ii=0; ii<pdfs.size( ); ii++)
    GetMostLikelyValue( CNval[ii], pdfs[ii] );
}


#endif
