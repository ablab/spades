///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Fastavector.h"
#include "paths/CAltFasta.h"
#include "paths/FixSomeIndelsUtils.h"
#include "efasta/EfastaTools.h"



/**
 * CAltFasta
 * Constructor
 */
CAltFasta::CAltFasta( int cid, int begin, int end, const fastavector &contig, 
    const vec<int>& sites,  const vec<vec<fastavector> >& sites_b, int wing1,
    int wing2 ) :
  contig_ ( &contig ),
  cid_ ( cid ),
  begin_ ( begin ),
  end_ ( end ),
  sites_ (sites),
  sites_b_ (sites_b),
  wing1_(wing1),
  wing2_(wing2),
  plog_ ( NULL )
{
  // sanity check
  ForceAssertEq( sites_.size(), sites_b_.size());
  // the original contig
  fastavector consensus;
  consensus.SetToSubOf( *contig_, begin, end - begin );
  // other alternatives
  vec<unsigned int> choice( sites_.size(), 0);
  do {
    // generate the new fastaString for that bases choices
    fastavector altFasta;
    for( size_t pos  = 0, k = 0; pos < consensus.size(); pos++ ) {
      if ( k < choice.size() && (int) pos == sites_[k] ) {
	altFasta = Cat( altFasta, sites_b_[k][choice[k]] );
	k++;
      } else altFasta.push_back( consensus[pos] );
    } 
    alt_.push_back( altFasta );
  } while( NextChoice(choice) );
  ForceAssert( consensus == alt_.front() );
}

///**
// * CAltFasta
// * AddAlternative
// */
//void CAltFasta::AddAlternative( const fastavector &alt )
//{
//  alt_.push_back( alt );
//}

/**
 * CAltFasta
 * SetOstringstream
 */
void CAltFasta::SetOstringstream( ostream &log )
{
  plog_ = &log;
}


/**
 * CAltFasta
 * EvalAlternatives
 *
 * A wrapper around several functions from David. The idea is to align
 * the reads belonging to the range [begin_, end_) on the contig cid_
 * on each possible alternative, and choose which alternatives are the
 * ones with better alignment scores.
 */
void CAltFasta::EvalAlternatives( const Bool VERBOSE,
				  const BaseVecVec &bases,
				  const BaseVecVec &jbases,
				  const VecQualNibbleVec &quals,
				  const VecQualNibbleVec &jquals,
				  const vec<tripal> &RALIGNS_this,
				  const vec<tripal> &JRALIGNS_this )
{
  // Heuristics.
  const int agree_size = 2;              // extra space on base sides
  // We allow as many error in the region. 
  const double min_core_identity = 0.0; 
  //const double min_core_identity = 1.0 - ( end_ - begin_ - sites_.size() ) * 0.10 / ( end_ - begin_ );  
  ForceAssert( reads_.empty() );

  //// Populate reads_ and readsq_.
  //SelectReads( VERBOSE, USE_QUALS, begin_, end_,
  //             RALIGNS_this, JRALIGNS_this, run_dir, TIGS,
  //             bases, jbases, quals, jquals, reads_, readsq_, plog_ );
  *plog_ << "quals.size()= " << quals.size() << endl;
  *plog_ << "jquals.size()= " << jquals.size() << endl;
  SelectReads2( VERBOSE, begin_, end_,
               RALIGNS_this, bases, quals, reads_, readsq_, plog_ );
  if ( VERBOSE )
    * plog_ << " Selected " << reads_.size() << " frag_reads" << endl;
  SelectReads2( VERBOSE, begin_, end_,
               JRALIGNS_this, jbases, jquals, reads_, readsq_, plog_ );
  if ( VERBOSE )
    * plog_ << " Selected " << reads_.size() << " total_reads" << endl;
 
  // Resize errs_.
  errs_.clear( );
  errs_.resize( reads_.size( ) );
  
  // Align reads to all alternatives (populate errs_).
  int cg_start = Max( 0, begin_ - wing1_ );
  int cg_end = Min( end_ + wing2_, (int)(*contig_).size( ) );
  for (int alt_id=0; alt_id<(int)alt_.size( ); alt_id++) {
    fvec left_chunk;
    fvec right_chunk;
    left_chunk.SetToSubOf( *contig_, cg_start, begin_ - cg_start );
    right_chunk.SetToSubOf( *contig_, end_, cg_end - end_ );
    fvec glued_chunks = Cat( left_chunk, alt_[alt_id], right_chunk );

    // I do not understand what agre_size is for.
    int nc1 = left_chunk.size( );
    int nc2  = alt_[alt_id].size();
    //int nc2  = right_chunk.size( );
    int NC1 = Max( 0, nc1 - agree_size );
    int NC2 = Min( nc2 + 2 * agree_size, (int)(*contig_).size( ) - NC1 );

    // Align reads.
    String title = "alt_" + ToString( alt_id );
    FindBestHomeForReads2( VERBOSE, title, reads_, readsq_,
			  glued_chunks, NC1, NC2, min_core_identity,
			  errs_, plog_ );
  }
 
}



void CAltFasta::Vote( bool VERBOSE )
{
  if (VERBOSE) *plog_ << "voting (alt:score)" << endl << endl;

  ostream *p_rout = plog_;
  vec< vec<String> > rows;
  for ( size_t i = 0; i < errs_.size(); i++ )
  {    
    vec<int> v = errs_[i];
    UniqueSort(v);
    if ( v.size( ) > 1 ) 
    {    
      const vec<int>& errVec = errs_[i];
      size_t n_err = errVec.size( );
      int best = Max(errVec);
      vec<double> SCORE(n_err);
      double total = 0.0;
      for ( size_t j = 0; j < n_err; j++ )
      {    
	if ( best - errVec[j] <= 10 )
        {    
	  SCORE[j] = pow( 10.0, double(errVec[j]-best)/10.0 );
          total += SCORE[j];    }    
      }
      for ( size_t j = 0; j < n_err; j++ )
      {    
	if ( best - errVec[j] <= 10 )
        {    
	  SCORE[j] /= total;
          votes.push( (signed)j, SCORE[j] );    
	}    
      }
      vec<String> row;
      if (VERBOSE) 
      {    
	row.push_back("[" + ToString(i) + "]");
        bool first = true;
        String scores, bests;
        for ( size_t j = 0; j < n_err; j++ ) 
        {    
	  if ( errVec[j] > -100000 )
          {    
	    if ( !first ) scores += ", "; 
            first = false;
            scores += ToString((signed)j) 
              + ":" + ToString(errVec[j]);    }
        }
        row.push_back(scores);
        first = True;
        for ( size_t j = 0; j < n_err; j++ ) 
        {    
	  if ( best - errVec[j] <= 10 )
          {    
	    if ( !first ) bests += ",";
            first = False;
            bests += ToString((signed)j) 
              + ":" + ToString(SCORE[j]);    }    
        }
        if ( bests.size( ) > 0 ) row.push_back(bests);
        rows.push_back(row);    
      }    
    }    
  }
  if (VERBOSE) PrintTabular(*p_rout, rows, 3);    
}

  
/// Decide which of the alternatives can be accepted based on the votes
void CAltFasta::SummarizeAndAcceptVotes( bool VERBOSE )
{
  ostream *p_rout = plog_;

  // prepare vote2 generated based on vote
  vec< pair<double,int> > votes2;
  Sort(votes);
  for (size_t i = 0; i < votes.size(); i++) 
  {    double X = 0.0;
    size_t j;
    for ( j = i; j < votes.size( ); j++ )
    { 
      if ( votes[j].first != votes[i].first ) break;
      X += votes[j].second;    
    }
    votes2.push( X, votes[i].first );
    i = j - 1;    
  }
  ReverseSort(votes2);
  if (VERBOSE) 
  {    *p_rout << "\nmax votes = " 
    << ( votes2.empty( ) ? 0 : votes2[0].first ) << "\n";
    vec< vec<String> > rows;
    vec<String> row;
    row.push_back( "votes", "alt" );
    rows.push_back(row);
    for (size_t i = 0; i < votes2.size(); i++) 
    {    vec<String> row;
      row.push_back(ToString(votes2[i].first), ToString(votes2[i].second));
      rows.push_back(row);    }
      *p_rout << "\n";
      PrintTabular(*p_rout, rows, 2);    
  }

  // Decide on change and announce.
  accepted.clear( );
  double vote_multiplier = 1.5;
  for (size_t i = 0; i < votes2.size(); i++) 
  {    
    if ( vote_multiplier * votes2[i].first < votes2[0].first ) break;
    accepted.push_back( votes2[i].second );    
  }    

  const double min_votes = 1.8;
  if ( accepted.solo( ) && accepted[0] == 0 ) return; // suggest the original sequence
  if ( accepted.empty( ) 
      || ( accepted.nonempty( ) && votes2[0].first < min_votes ) )
  {  
    accepted.clear( );
    accepted.push_back( 0 );    
  }
  if (VERBOSE) 
  {    *p_rout << "\nRECOMMEND CHANGE:";
    for (int i = 0; i < accepted.isize(); i++)
      *p_rout << " " << accepted[i];
    *p_rout << "\n";    
  }
  Sort(accepted);

  // If the form present in the original contig is in accepted,
  // move it to the front.
  if ( Member( accepted, 0 ) )
  {    
    vec<int> accepted2;
    accepted2.push_back(0);
    for ( int l = 0; l < accepted.isize( ); l++ )
    {    if ( accepted[l] != 0 ) 
      accepted2.push_back( accepted[l] );    
    }
    accepted = accepted2;    
  }
}


/// Generate all accepted edits, which is a vector of triples of (start, stop, replacement)
void CAltFasta::AcceptedEdits(vec< triple<int, int, String> > &edits, bool VERBOSE)
{
  if ( accepted.empty() ) return;
  if ( accepted.solo() && accepted.front() == 0 ) return;
  // which alternatives has been taken?
  vec<vec<bool> > flagPick(sites_.size() );
  for( size_t i = 0; i < sites_.size(); i++ )
    flagPick[i].resize( sites_b_[i].size(), false);
  // now mark them
  {
    vec<unsigned int> choice( sites_.size(), 0 );
    int index = 0; 
    do {
      if ( find( accepted.begin(), accepted.end(), index ) != accepted.end() )
	for ( size_t k = 0; k < choice.size(); k++ )
	  flagPick[k][choice[k]] = true;
      index++;
    }
    while ( NextChoice(choice) );
  }
  // generate the efasta string
  const fastavector& original = alt_[0];
  for( size_t k = 0; k < sites_.size(); k++ )
  {
    set<fastavector> temp;
    for ( size_t kk = 0; kk < flagPick[k].size(); kk++)
      if ( flagPick[k][kk] ) 
	temp.insert( sites_b_[k][kk] );
    if ( temp.size() == 0 ) continue;
    // prepare the efasta string
    EFastaVec altFasta;
    if (temp.size() == 1 )
    {
      if ( *temp.begin() ==  sites_b_[k][0] ) continue; // no change
      if ( temp.begin()->size() > 0 ) 
	altFasta.append( temp.begin()->ToString() );
    }
    else
    {
      altFasta.push_back('{');
      unsigned int counter = 0;
      for( set<fastavector>::iterator it = temp.begin(); it != temp.end(); ++it )
      {
	if ( it->size() > 0 )
	  altFasta.append( it->ToString() );
	if ( ++counter < temp.size()) altFasta.push_back(',');
      }
      altFasta.push_back('}');
      // sy_fix_later
      // only if unambiguously determined
      if ( VERBOSE ) *plog_ << "  Ambiguous Fasta skipped" << altFasta << endl;
      continue;
    }
    // output
    edits.push( begin_ + sites_[k], begin_ + sites_[k] + 1 , altFasta );
    if (VERBOSE)
    {
      *plog_ << "OldFasta " << original[sites_[k]] << endl;
      *plog_ << "NewFasta " << altFasta << endl;
    }
  } // end for each site

}

// Print the debug information
void CAltFasta::Print( ostream & out )
{
  out << "\n======== CAltFasta debug: " << cid_ << " " << begin_ << " " << end_ << " At sites: " ;
  for ( size_t i = 0; i < sites_.size(); i++ ) out << sites_[i] << " ";
  out << " wings: " << wing1_ << "," << wing2_ << "  =======" << endl;
  out << "Alternatives " << alt_.size() << endl;
  String altSites( alt_[0].size(), ' ');
  //for ( size_t k = 0; k < sites_.size(); k++)
  //  altSites[sites_[k]] = '*';
  for( size_t i = 0; i < alt_.size(); i++ )
  {
    if ( ! alt_[i].size() == 0  )
      out << i << " " << alt_[i].ToString() << endl;
    else
      out << i << "  "  << endl;
  }
  out << "  " << altSites << endl;
  out << "reads_.size()= " << reads_.size() << endl;
  out << "votes.size()= " << votes.size() << endl;
  //for( size_t i = 0; i < votes.size(); i++ )
  //  out << "votes " << votes[i].first << ":" << votes[i].second << endl;
  for( size_t i = 0; i < accepted.size(); i++ )
    out << "accepted " << accepted[i] << endl;
  //for( size_t i = 0; i < errs_.size(); i++ )
  //{
  //  out << "errs " << i << " " << errs_[i].size() << ":";
  //  for( size_t j = 0; j < errs_[i].size(); j++ )
  //    out << errs_[i][j] << " ";
  //  out << endl;
  //}
  out << endl;
}


// Return the next enumeration of a vector "choice" where sites_[i] will
// pick the alternative fasta sites_b_[i][choice[i]]
bool CAltFasta::NextChoice(vec<unsigned int>& choice)
{                                                                                                                                                                                                                                                                       
  int loc = choice.size() -1;                                                                                                                                                                                                                                           
  while(loc >= 0)                                                                                                                                                                                                                                                       
    if ( choice[loc] < sites_b_[loc].size() - 1 )                                                                                                                                                                                                                                   
    {                                                                                                                                                                                                                                                                   
      choice[loc]++;                                                                                                                                                                                                                                                    
      return true;                                                                                                                                                                                                                                                      
    }                                                                                                                                                                                                                                                                   
    else{                                                                                                                                                                                                                                                               
      choice[loc] = 0;                                                                                                                                                                                                                                                  
      loc--;                                                                                                                                                                                                                                                            
    }                                                                                                                                                                                                                                                                   
  return false;                                                                                                                                                                                                                                                         
}              
