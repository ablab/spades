///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PATHS__C_ALT_FASTA
#define PATHS__C_ALT_FASTA

#include "Basevector.h"
#include "Fastavector.h"
#include "feudal/QualNibbleVec.h"

typedef triple<int64_t,int,int> tripal;

/**
 * class CAltFasta
 *
 * It contains multiple possible variants of a given interval on a
 * genomic region. This class is used to locally fix consensus, and it
 * encapsulates David's code originally in FixSomeIndels.
 *
 * How it works: instantiate a CAltFasta object, and then add as many
 * alternate forms as you need. EvalAlternatives will return a set of
 * fixes.
 */
class CAltFasta {
  
public:
  
  CAltFasta( int cid, int begin, int end, const fastavector &contigs, 
      const vec<int>& sites, const vec<vec<fastavector> >& sites_b, 
      int wing1 = 100, int wing2 = 100 );
  
  // Add an alternative sequence to consensus.
  //void AddAlternative( const fastavector &alternative );

  // Thread-safe log stringstream.
  void SetOstringstream( ostream &log );
  
  // Evaluate alternatives, return corrections.
  void EvalAlternatives( const Bool VERBOSE,
			 const BaseVecVec &bases,
			 const BaseVecVec &jbases,
			 const VecQualNibbleVec &quals,
			 const VecQualNibbleVec &jquals,
			 const vec<tripal> &RALIGNS_this,
			 const vec<tripal> &JRALIGNS_this );
  
  // Decide the number of votes for each alternative sequence
  void Vote( bool VERBOSE = true );
  
  /// Decide which of the alternatives can be accepted based on the votes
  void SummarizeAndAcceptVotes( bool VERBOSE = False );
  /// Generate all accepted edits, which is a vector of triples of (start, stop, replacement)
  void AcceptedEdits(vec< triple<int, int, String> > &edits, bool VERBOSE = False);
  
  // Print the debug information
  void Print( ostream & out );
  
  
private:
  
  const fvec *contig_;  // fastavector of contig
  int cid_;             // id of contig
  int begin_;           // begin of interval on contig
  int end_;             // end of interval

  vec<int> sites_; // identify the locations in the original contigs that needs to be changed
                   // begin_ + site[i] is the exact position in the contig
  vec<vec<fastavector> > sites_b_; // alternative bases at each sites, sites_b_[i] lists all possible
                                   // fastavectors of sites[i]. sites_b_[i][0] <- contig_[begin_ + sites[i]]

  int wing1_, wing2_;   // safe range on both sides

  vecfvec alt_;           // alternates (alt_[0] = [begin, end) in contig_)
  vecbvec reads_;         // reads over [begin, end) of contig cid_
  vecqnibvec readsq_;     // qualities of reads
  vec< vec<int> > errs_;  // errs_[i][j]: score of align reads_[i] on alt_[j]
  ostream *plog_;   // thread-safe log stream

  vec< pair<int,double> > votes;
  vec<int> accepted;
  
private: // other private methods
  // Return the next enumeration of a vector "choice" where sites_[i] will
  // pick the alternative fasta sites_b_[i][choice[i]]
  bool NextChoice(vec<unsigned int>& choice);

};

#endif
