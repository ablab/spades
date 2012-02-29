///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#include "AnnotatedContig.h"
#include "system/System.h"

const vec< basevector >*          semiannotation::standards_         = NULL;
const vec< basevector >*          annotation::standards_             = NULL;
const vec< basevector >*          annotated_contig::standards_       = NULL;
const vec< basevector >*          annotated_supercontig::standards_  = NULL;
const map< int, arachne_contig >* annotated_contig::arachne_contigs_ = NULL;

bool semiannotation_loccomp::operator() ( const semiannotation &ca1,
					  const semiannotation &ca2 ) const {

  if      ( ca1.RC() != ca2.RC() )
    return ( ca2.RC() );
  else if ( ( ca1.StartOnStandard() - ca1.StartOnThis() ) !=
	    ( ca2.StartOnStandard() - ca2.StartOnThis() ) )
    return ( ( ca1.StartOnStandard() - ca1.StartOnThis() ) <
	     ( ca2.StartOnStandard() - ca2.StartOnThis() ) );
  else if ( ca1.StartOnThis() != ca2.StartOnThis() )
    return ( ca1.StartOnThis() < ca2.StartOnThis() );
  else
    return false;
}


int arachne_contig_placement::ID() const {
  if ( arachne_id_.Contains( "." ) ) {

    String w = arachne_id_.Before( "." );

    return ( atoi( w.c_str() ) );
  }
  else
    return ( atoi( arachne_id_.c_str() ) );
}

void arachne_contig_placement::Print( ostream &o ) const {
  o << "ID: " << ID() << "; [from,to]: ["
    << start_ << ","
    << stop_ << "]";
  if ( IsGap() )
    o << " -- GAP";
}

void annotated_contig::ComputeSemiannotations() {
  if ( ! NumACPs() ) {
    SetSemiannotated( False );
    return;
  }
  acp_itr curracp = FirstACP();
  acp_itr end_acp = EndACP();

  typedef set< semiannotation, semiannotation_loccomp > set_loccomp;
  typedef set_loccomp::iterator                         set_loccomp_itr;
  typedef map< int, set_loccomp >::iterator             map_set_loccomp_itr;

  map< int, set_loccomp > stdIDs_and_locs; /* The 'key' of this map is the standard id,
					      and the set< contig_actualloc, ctg_actloc_loccomp >
					      contig_actuallocs, sorted according to location,
					      found for
					      this annotated_contig, by looking at all the different
					      arachne contig actuallocs for the non-gap
					      arachne_contig_placements in this annotated_contig. */

  // Print( cout );
  // cout << endl;
  for ( ;
	curracp != end_acp;
	++curracp ) {
    if ( curracp->IsGap() )
      continue;

    const arachne_contig &arach_ctg = ArachneContig( curracp->ID() );

    // cout << "\t\tArachne Contig  ";
    // arach_ctg.Print( cout );
    // cout << endl;

    ctg_actloc_itr curractloc = arach_ctg.FirstActloc();
    ctg_actloc_itr end_actloc = arach_ctg.EndActloc();

    for ( ;
	  curractloc != end_actloc;
	  ++curractloc ) {

      // cout << "\t\t\t Actloc: ";
      // curractloc->Print( cout );
      // cout << endl;

      int implied_left_on_std = ( curractloc->RC() ?
				  curractloc->Start() - Length() + curracp->Start() + arach_ctg.Length() :
				  curractloc->Start() - curracp->Start() );

      int second_implied_left_on_std = ( curractloc->RC() ?
					 curractloc->Start() - Length() + curracp->Stop() :
					 curractloc->Start() + arach_ctg.Length() - curracp->Stop() );

      for ( int i_choose = 0;
	    i_choose < 2;
	    ++i_choose ) {

	if ( i_choose )
	  if ( abs( implied_left_on_std - second_implied_left_on_std ) < max_fudge / 2 )
	    continue;

	int chosen_implied_left_on_std = ( i_choose
					   ? second_implied_left_on_std
					   : implied_left_on_std );

	int chosen_left_on_std = ( i_choose
				   ? curractloc->Start() + second_implied_left_on_std - implied_left_on_std
				   : curractloc->Start() );


	semiannotation this_semi( curractloc->StdID  (),
				  (unsigned char)curractloc->RC(),
				  curracp->Start     (),
				  chosen_left_on_std,
				  chosen_implied_left_on_std,
				  arach_ctg.Length   (),
				  curractloc->Weight () );

        /*
	cout << "\t\t\t Semi: ";
	this_semi.Print( cout );
	cout << endl
	     << "\t\t\t\t Implied boundaries on standard: ( " << chosen_implied_left_on_std
	     << " -> " << chosen_implied_left_on_std + Length()
	     << " )" << " two options: " << implied_left_on_std << "," << second_implied_left_on_std
	     << endl;
        */

	map_set_loccomp_itr curr_stdID_set_loccomp_itr = stdIDs_and_locs.find( curractloc->StdID() );

	if ( curr_stdID_set_loccomp_itr != stdIDs_and_locs.end() ) {

	  Assert( curr_stdID_set_loccomp_itr->first == curractloc->StdID() );

	  set_loccomp_itr find_this_semi = curr_stdID_set_loccomp_itr->second.find( this_semi );

	  if ( find_this_semi != curr_stdID_set_loccomp_itr->second.end() ) {
	    this_semi.AddToWeight( find_this_semi->Weight() );
	    // cout << "\t\t\t\t ";
	    // this_semi.Print( cout );
	    // cout << " == ";
	    // find_this_semi->Print( cout );
	    // cout << endl;

	    curr_stdID_set_loccomp_itr->second.erase( find_this_semi );
	    curr_stdID_set_loccomp_itr->second.insert( this_semi );

	    // cout << "\t\t\t\t Added to the weight of this semi: ";
	    // this_semi.Print( cout );
	    // cout << endl;

	    set_loccomp_itr currsemi_print = curr_stdID_set_loccomp_itr->second.begin();
	    set_loccomp_itr end_semi_print = curr_stdID_set_loccomp_itr->second.end();

	    for ( ;
		  currsemi_print != end_semi_print;
		  ++currsemi_print ) {
	      // cout << "\t\t\t\t\t";
	      // currsemi_print->Print( cout );
	      // cout << endl;
	    }
	  }
	  else {
	    curr_stdID_set_loccomp_itr->second.insert( this_semi );
	    // cout << "\t\t\t\t Inserted this semi" << endl;
	  }
	}
	else {
	  stdIDs_and_locs.insert( pair< int, set_loccomp > ( curractloc->StdID(),
							     set_loccomp() ) );
	  curr_stdID_set_loccomp_itr = stdIDs_and_locs.find( curractloc->StdID() );

	  Assert( curr_stdID_set_loccomp_itr        != stdIDs_and_locs.end() );
	  Assert( curr_stdID_set_loccomp_itr->first == curractloc->StdID()   );

	  curr_stdID_set_loccomp_itr->second.insert( this_semi );
	  // cout << "\t\t\t\t Inserted this semi, after inserting STANDARD" << endl;
	}

	// cout << "\t\t\t\t\t Semi count: " << curr_stdID_set_loccomp_itr->second.size() << endl;
      }
    }
  }

  map_set_loccomp_itr currset_of_locs_itr = stdIDs_and_locs.begin();
  map_set_loccomp_itr end_set_of_locs_itr = stdIDs_and_locs.end();

  for ( ;
	currset_of_locs_itr != end_set_of_locs_itr;
	++currset_of_locs_itr ) {

    int          curr_std_id     = currset_of_locs_itr->first;
    set_loccomp &currset_of_locs = currset_of_locs_itr->second;

    int fudge = MIN( max_fudge,
		     ( (int)bases_.size() / fudge_per_nucleotides ) + constant_fudge );

    set_loccomp_itr currsemi_itr = currset_of_locs.begin();
    set_loccomp_itr end_semi_itr = currset_of_locs.end();

    for ( ;
	  currsemi_itr != end_semi_itr;
	  /* ++currsemi_itr */ ) {
      /*
      cout << "\t curr semi: ";
      currsemi_itr->Print( cout );
      cout << endl
           << "\t\t\t\t Implied boundaries on standard: ( " << currsemi_itr->ShiftFromStandard()
	   << " -> " << currsemi_itr->ShiftFromStandard() + Length()
	   << " )" << endl;
      */
      if ( fudge + currsemi_itr->StartOnStandard() < 0 ) {
	// cout << "\t\t\t\t BAD SEMI!" << endl;
	++currsemi_itr;
	continue;
      }

      set_loccomp_itr prevsemi_itr = currsemi_itr;
      longlong total_actlocs  = prevsemi_itr->Weight();
      longlong total_shift    = (longlong)prevsemi_itr->Weight() * (longlong)prevsemi_itr->ShiftFromStandard();
      longlong minstart       = prevsemi_itr->StartOnThis();
      longlong maxstop        = prevsemi_itr->StartOnThis() + prevsemi_itr->Length() - 1;

      /*
      cout << "\t\t\t\t tot_actloc: " << total_actlocs
	   << "; tot_shift: " << total_shift
	   << ", " << total_shift/ total_actlocs + minstart
	   << "; minstart,maxstop: " << minstart << ", " << maxstop << endl;
      */
      for ( ++currsemi_itr;
	    currsemi_itr != end_semi_itr &&
	      currsemi_itr->RC() == prevsemi_itr->RC() &&
	      ( currsemi_itr->ShiftFromStandard() - prevsemi_itr->ShiftFromStandard() < fudge);
	    ++currsemi_itr ) {

        /*
  	cout << "\t Curr semi: ";
	currsemi_itr->Print( cout );
	cout << endl
	     << "\t\t\t\t Implied boundaries on standard: ( " << currsemi_itr->ShiftFromStandard()
	     << " -> " << currsemi_itr->ShiftFromStandard() + Length()
	     << " )" << endl;
        */

	++prevsemi_itr;

	Assert( prevsemi_itr == currsemi_itr );

	total_actlocs += prevsemi_itr->Weight();
	total_shift   += (longlong)prevsemi_itr->Weight() * (longlong)prevsemi_itr->ShiftFromStandard();
	minstart       = MIN( prevsemi_itr->StartOnThis(), minstart );
	maxstop        = MAX( prevsemi_itr->StartOnThis() + prevsemi_itr->Length() - 1, maxstop );

        /*
	cout << "\t\t\t\t tot_actloc: " << total_actlocs
	     << "; tot_shift: " << total_shift
	     << ", " << total_shift/ total_actlocs + minstart
	     << "; minstart,maxstop: " << minstart << ", " << maxstop << endl;
        */
      }

      /*
      cout << "\t\t\t\t\t\t Exited loop because ";
      if ( currsemi_itr == end_semi_itr ) cout << " reached end of semiannotations " << endl;
      else if ( currsemi_itr->RC() != prevsemi_itr->RC() ) cout << " conflicting RC bits" << endl;
      else cout << " too little fudge " << endl;

      cout << "\t\t\t\t tot_actloc: " << total_actlocs
	   << "; tot_shift: " << total_shift
	   << ", " << total_shift/ total_actlocs + minstart
	   << "; minstart,maxstop: " << minstart << ", " << maxstop << endl;
      */


      longlong implied_length      = min( (longlong)Length(), maxstop - minstart + 1 );
      longlong implied_minstart    = ( prevsemi_itr->RC() ?
					    max( (longlong)0, Length() - maxstop ):
					    minstart );

      longlong start_on_std        = total_shift / total_actlocs + implied_minstart;
      longlong implied_left_on_std = start_on_std - implied_minstart;


      /*
      cout << "\t\t\t\t Implied left_on_std, length, minstart, maxstop: " << implied_left_on_std << ", " << implied_length << ", ( "
	   << implied_minstart << " , " << implied_minstart + implied_length - 1 << " )" << endl;
      */

      semiannotation new_semi( curr_std_id,
			       prevsemi_itr->RC(),
			       implied_minstart,
			       start_on_std,
			       implied_left_on_std,
			       implied_length,
			       total_actlocs );

      // cout << "\t  NEW semi: ";
      // new_semi.Print( cout );
      // cout << endl;
      semiannotations_.push_back( new_semi );
      if ( currsemi_itr == end_semi_itr )
	break;
    }
  }
  SetSemiannotated( True );
}

bool annot_sc_lengthcomp ( const annotated_supercontig &a,
			   const annotated_supercontig &b ) {
  return a.Length() > b.Length();
}

int annotated_supercontig::Length() const {
  int length = 0;

  for ( int i = 0;
	i < (int)contigs_.size();
	++i )
    length += contigs_[ i ].Length();

  return length;
}

void annotated_supercontig::ComputeSemiannotations() {
  annotated_contig_itr currctg_itr = contigs_.begin();
  annotated_contig_itr end_ctg_itr = contigs_.end();

  for ( ;
	currctg_itr != end_ctg_itr;
	++currctg_itr )
    currctg_itr->ComputeSemiannotations();
}

void annotated_final_answer::ComputeSemiannotations() {
  annotated_supercontig_itr currsc_itr = supers_.begin();
  annotated_supercontig_itr end_sc_itr = supers_.end();

  for ( ;
	currsc_itr != end_sc_itr;
	++currsc_itr )
    currsc_itr->ComputeSemiannotations();
}


ostream& operator<<( ostream &o, const semiannotation &a ) {
  o << a.std_id_                  << " "
    << BoolToInt( a.rc_ )         << " "
    << a.start_on_this_           << " "
    << a.start_on_std_            << " "
    << a.implied_leftmost_on_std_ << " "
    << a.length_                  << " "
    << a.weight_                  << " ";

  return o;
}

void semiannotation::Print( ostream &o ) const {
  o << "start,stop: ("      << StartOnThis()
    << ","       << (StartOnThis() + Length() - 1)
    << ")  -->  stdID: " << StdID()
    << "; start,stop,dir,shift on std.: (" << StartOnStandard()
    << ","       << ( StartOnStandard() + Length() - 1 )
    << "), "    << ( RC() ? "rc" : "fw" )
    << " " << ShiftFromStandard()
    << ";  wht: "      << Weight();
}

void annotated_contig::Print( ostream &o ) const {
  o << "Contig " << id_ << ": length = " << quals_.size() << ".";
}

void annotated_supercontig::Print( ostream &o ) const {
  o << "Supercontig of length " << Length() << ":" << endl;
  for ( int i = 0; i < (int)contigs_.size(); ++i ) {
    o << "\t";
    contigs_[ i ].Print( o );
    o << endl;
  }
}



ostream& operator<<( ostream &o, const annotation &a ) {
  o << a.std_id_          << " "
    << BoolToInt( a.rc_ ) << " ";
  a.align_.Write( o );
  o << " ";
  return o;
}

ostream& operator<<( ostream &o, const arachne_contig_placement &a ) {
  o << a.arachne_id_ << " "
    << a.start_      << " "
    << a.stop_       << " ";

  return o;
}


ostream& operator<<( ostream &o, const annotated_contig &a ) {
  AssertEq( a.bases_.size(), a.quals_.size() );

  o << a.ID() << " ";
  o << a.bases_.size() << " ";
  for ( int i = 0; i < (int)a.bases_.size(); ++i )
    o.put( a.bases_[ i ] );

  for ( int i = 0; i < (int)a.quals_.size(); ++i )
    o.put( a.quals_[ i ] );

  o << '\n';

  o << a.semiannotations_.size() << endl;
  for ( int i = 0; i < (int)a.semiannotations_.size(); ++i )
    o << a.semiannotations_[ i ] << endl;

  o << a.annotations_.size() << endl;
  for ( int i = 0; i < (int)a.annotations_.size(); ++i )
    o << a.annotations_[ i ] << endl;

  o << a.arachne_ctg_placements_.size() << endl;
  for ( int i = 0; i < (int)a.arachne_ctg_placements_.size(); ++i )
    o << a.arachne_ctg_placements_[ i ] << endl;

  return o;
}




ostream& operator<<( ostream &o, const annotated_supercontig &a ) {
  Assert( a.contigs_.size() == a.gaps_.size() + 1 &&
	  a.gaps_.size()    == a.gap_sds_.size() );


  o << a.contigs_.size() << endl;
  for ( int i = 0; i < (int)a.contigs_.size(); ++i )
    o << a.contigs_[ i ] << endl;

  o << a.gaps_.size() << endl;

  for ( int i = 0; i < (int)a.gaps_.size(); ++i )
    o << a.gaps_[ i ] << " " << a.gap_sds_[ i ] << " ";

  return o;
}



ostream& operator<<( ostream &o, const annotated_final_answer &a ) {

  o << a.supers_.size() << endl;

  for ( int i = 0; i < (int)a.supers_.size(); ++i )
    o << a.supers_[ i ] << endl;

  return o;
}


istream& operator>>( istream &i, semiannotation &a ) {
  i >> a.std_id_;
  int dummy;
  i >> dummy;
  a.rc_ = ( Bool ) dummy;

  i >> a.start_on_this_;
  i >> a.start_on_std_;
  i >> a.implied_leftmost_on_std_;
  i >> a.length_;
  i >> a.weight_;


  Assert( a.std_id_        >= 0 );
  Assert( a.start_on_this_ >= 0 );
  Assert( a.start_on_std_  >= 0 );
  Assert( a.length_        >= 0 );
  Assert( a.weight_        >= 0 );

  return i;
}


istream& operator>>( istream &i, annotation &a ) {
  i >> a.std_id_;
  int dummy;
  i >> dummy;
  a.rc_ = ( Bool ) dummy;
  a.align_.Read( i );

  Assert( a.std_id_ >= 0 );
  return i;
}

istream& operator>>( istream &i, arachne_contig_placement &a ) {
  i >> a.arachne_id_
    >> a.start_
    >> a.stop_;

  /*
  if ( !( a.start_  >= 0 && a.stop_ > a.start_ ) )
  {    cout << "\nillegal arachne contig placement detected\n";
       PRINT3( a.arachne_id_, a.start_, a.stop_ );
       exit(1);     }
  */

  return i;
}

istream& operator>>( istream &in, annotated_contig &a ) {
  int id;
  in >> id;
  a.SetID( id );

  int len;
  in >> len;
  Assert( len >= 0 );

  a.bases_.Setsize( len );
  a.quals_.resize( len );

  char dummy_char;
  in.get( dummy_char );
  Assert( dummy_char == ' ' );
  for ( int i = 0; i < len; ++i ) {
    in.get( dummy_char );
    a.bases_.Set( i, (unsigned char) dummy_char );
  }

  for ( int i = 0; i < len; ++i )
  {
    in.get( dummy_char );
    a.quals_[ i ] = dummy_char;
  }

  in.get( dummy_char );
  Assert( dummy_char == '\n' );
  in.putback( dummy_char );

  in >> len; // a.semiannotations_.size()
  Assert( len >= 0 );
  a.semiannotations_.resize( len );

  for ( int i = 0; i < len; ++i )
    in >> a.semiannotations_[ i ];

  in >> len; // a.annotations_.size()
  Assert( len >= 0 );
  a.annotations_.resize( len );

  for ( int i = 0; i < len; ++i )
    in >> a.annotations_[ i ];

  in >> len; // a.annotations_.size()
  Assert( len >= 0 );
  a.arachne_ctg_placements_.resize( len );
  for ( int i = 0; i < len; ++i )
    in >> a.arachne_ctg_placements_[ i ];


  return in;
}

istream& operator>>( istream &in, annotated_supercontig &a ) {
  int len;

  in >> len;
  Assert( len >= 0 );
  a.contigs_.resize( len );

  for ( int i = 0; i < len; ++i )
    in >> a.contigs_[ i ];

  in >> len;
  Assert( len >= 0 );
  a.gaps_.resize   ( len );
  a.gap_sds_.resize( len );

  Assert( a.gaps_.size() + 1 == a.contigs_.size() );
  for ( int i = 0; i < len; ++i )
    in >> a.gaps_[ i ] >> a.gap_sds_[ i ];

  return in;
}


istream& operator>>( istream &in, annotated_final_answer &a ) {

  int len;
  in >> len;
  Assert( len >= 0 );

  for ( int i = 0; i < len; ++i ) {
    annotated_supercontig asc;
    in >> asc;
    a.PushSupercontig( asc );
  }

  return in;
}

arachne_contig& annotated_final_answer::ArachneContig( int i ) {
  arachne_ctg_map_itr acm_itr = arachne_contigs_.find( i );
  if ( acm_itr == arachne_contigs_.end() )
    crash();
  return acm_itr->second;
}

