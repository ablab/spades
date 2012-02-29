// Copyright (c) 2000, 2001 Whitehead Institute for Biomedical Research
// 

#include "agp/AgpFile.h"
#include "TokenizeString.h"

char agp_contig::TypeChar_() const
{
    switch ( type_ )
    {
        case wgs_contig:
            return 'W';
        case finished_contig:
            return 'F';
        case draft_contig:
            return 'D';
        case blunt_joined_contig:
            return 'A';
        default:
            return '0';
    }
}

agp_contig::contig_type TypeFromChar( const char c )
{
    switch ( c )
    {
        case 'W':
            return agp_contig::wgs_contig;
        case 'F':
            return agp_contig::finished_contig;
        case 'D':
            return agp_contig::draft_contig;
        case 'N':
            return agp_contig::blunt_joined_contig;
        default:
            return agp_contig::other_contig;
    }
}

bool agp_contig::Check( int cg_length ) const
{
  if ( this->Start( ) < 0 )
    return false;
  
  if ( this->Stop( ) > cg_length - 1 )
    return false;

  return true;
}

void agp_contig::Print( ostream& out,
			const String& chromosome,
			longlong& base_position,
			longlong& unit_position,
                        const String *base_name )
{
  int length = this->Extent();
  out << ( base_name ? *base_name : "chr" ) << chromosome << "\t"
      << base_position << "\t" 
      << base_position + length - 1 << "\t" 
      << unit_position << "\t" 
      << this->TypeChar_() << "\t"
      << this->Name() << "\t"  
      << this->Start() + 1 << "\t" 
      << this->Stop() + 1<< "\t";
  if ( this->IsRC() )
    out << "-";
  else
    out << "+";
  out << "\n";

  base_position += length;
  ++unit_position;
}

void agp_gap::Print( ostream& out,
		     const String& chromosome,
		     longlong& base_position,
		     longlong& unit_position,
                     const String *base_name )
{
  out << ( base_name ? *base_name : "chr" ) << chromosome << "\t"
      << base_position << "\t" 
      << base_position + this->Length() - 1 << "\t" 
      << unit_position << "\t" 
      << "N" << "\t"
      << this->Length() << "\t" 
      << this->Type() << "\t" 
      << ( this->IsBridged() ? "yes" : "no" ) << "\t\n";

  base_position += this->Length();
  ++unit_position;
}


/*
 * GetEntry
 *
 * The ii entry_type can be either a contig or a gap. If it is a contig,
 * then *contig will be set to the pointer pointing at a agp_contig, and
 * *gap to null. If instead ii is a gap, then *gap will be set to the
 * pointer pointing at the proper gap, and *contig to null.
 */
void agp_chromosome::GetEntry( int ii, agp_contig **contig, agp_gap **gap )
{
  entry_type type = entries_[ii].type;
  int index = entries_[ii].index;
  
  if ( contig_type == type ) {
    *gap = 0;
    *contig = &( contigs_[index] );
  }
  else {
    *contig = 0;
    *gap = &( gaps_[index] );
  }

}



// agp_chromosome
// 
// the entries_ member contains pairs binding the type of entry to the
// index in the appropriate vector.

void agp_chromosome::AddContig( const agp_contig& contig )
{
  contigs_.push_back( contig );
  entries_.push_back( entry_pair( contig_type, contigs_.size()-1 ) );
}

void agp_chromosome::AddGap( const agp_gap& new_gap, bool strict )
{
    // If the last entry was a gap...
    if ( ! entries_.empty() &&
         entries_.back().type == gap_type )
    {
        agp_gap& existing_gap = gaps_[entries_.back().index];

        existing_gap.SetLength( existing_gap.Length() + new_gap.Length() );

        // If the new gap is a clone gap and the existing gap is a
        // fragment gap, change the existing gap to a clone gap.
        if ( existing_gap.IsBridged() &&
             ! new_gap.IsBridged() )
        {
            existing_gap.SetType( "clone" );
            existing_gap.SetBridged( false );
        }
    }            
    else
    {
        gaps_.push_back( new_gap );
        entries_.push_back( entry_pair( gap_type, gaps_.size()-1 ) );
    }
}
  
void agp_chromosome::AddSuper( const vec<int>& contig_ids, 
			       const vec<int>& contig_lens, 
			       const vec<int>& gaps,
			       const int minimum_gap,
			       const agp_mods& modifications,
			       const bool is_reversed )
{
  AssertEq( contig_ids.size(), contig_lens.size() );
  AssertEq( gaps.size(), contig_ids.size() - 1 );

  // Local copies
  vec<int> ids_l; ids_l.reserve( contig_ids.size() );
  vec<int> lens_l; lens_l.reserve( contig_lens.size() );
  vec<int> gaps_l; gaps_l.reserve( gaps.size() );

  // Fill out the local copies, removing unplaced contigs and
  // adjusting or removing the preceding gaps if necessary.
  const vec<int>& contigs_to_unplace = modifications.GetContigsToUnplace();

  for ( unsigned int ii = 0; ii < contig_ids.size(); ++ii )
  {
    int contig_id = contig_ids[ii];

    // should we unplace this contig?
    if ( binary_search( contigs_to_unplace.begin(),
			contigs_to_unplace.end(),
			contig_id ) )
    {
      if ( ! gaps_l.empty() )
      {
	// if this isn't the last contig, then add the unplaced
	// contig's length and the length of the following gap to the
	// preceding gap.
	if ( ii < contig_ids.size() - 1)
	{
	  gaps_l.back() += contig_lens[ii];
	  gaps_l.back() += gaps[ii];
	}

	// if this is the last contig, remove the preceding gap.
	else
	  gaps_l.resize( gaps_l.size() - 1 );
      }
    }
    else
    {
      ids_l.push_back( contig_id );
      lens_l.push_back( contig_lens[ii] );
      if ( ii < gaps.size() )
	gaps_l.push_back( gaps[ii] );
    }
  }

  // Trim contigs, adjusting gaps if necessary.
  const vec<contig_trim>& trims = modifications.GetTrims();

  vec<int> starts( lens_l.size() );
  vec<int> stops( lens_l.size() );

  for ( unsigned int ii = 0; ii < ids_l.size(); ++ii )
  {
    starts[ii] = 0;
    stops[ii] = lens_l[ii] - 1;

    int contig_id = ids_l[ii];

    for ( unsigned int jj = 0; jj < trims.size(); ++jj )
    {
      if ( trims[jj].Id() != contig_id )
	continue;

      if ( trims[jj].TrimFront() )
      {
	starts[ii] += trims[jj].Amount();
	
	// adjust preceding gap
	if ( ii > 1 && ii <= gaps_l.size() )
	  gaps_l[ii-1] += trims[jj].Amount();
      }

      else
      {
	stops[ii] -= trims[jj].Amount();
	
	// adjust following gap
	if ( ii < gaps_l.size() )
	  gaps_l[ii] += trims[jj].Amount();
      }
    }
  }

  // Reverse the order of the supercontig, if necessary.
  if ( is_reversed )
  {
    ids_l.ReverseMe();
    lens_l.ReverseMe();
    starts.ReverseMe();
    stops.ReverseMe();
    gaps_l.ReverseMe();
  }

  AssertEq( ids_l.size(), lens_l.size() );
  AssertEq( gaps_l.size(), ids_l.size() - 1 );

  for ( unsigned int ii = 0; ii < ids_l.size(); ++ii )
  {    
    if ( ii > 0 )
    {
      int reported_gap = max<int>( gaps_l[ii-1], minimum_gap );
      this->AddGap( agp_gap( "fragment", reported_gap, true ) );
    }

    AssertLt( starts[ii], stops[ii] );

    this->AddContig( agp_contig( ids_l[ii], 
				 lens_l[ii],
				 starts[ii], stops[ii], 
				 is_reversed ) );
  }
}



bool agp_chromosome::AttemptInsertion( agp_contig& new_contig,
				       insertion ins,
				       const int minimum_gap_size )
{
  int first_anchor_entry_idx = -1;
  for ( unsigned int ii = 0; ii < entries_.size(); ++ii )
  {
    if ( entries_[ii].type == contig_type )
    {
      int contig_idx = entries_[ii].index;
      agp_contig& ctg = contigs_[ contig_idx ];

      if ( ins.IsFirstAnchor( ctg.Id() ) )
      {
        first_anchor_entry_idx = ii;
        break;
      }

      else if ( ins.IsLastAnchor( ctg.Id() ) )
      {
        first_anchor_entry_idx = ii;
        new_contig.SetRC( true );
        ins.Reverse();
        break;
      }
    }
  }

  if ( first_anchor_entry_idx == -1 )
    return false;

  int first_anchor_contig_idx = entries_[ first_anchor_entry_idx ].index;
  agp_contig& first_anchor_contig = contigs_[first_anchor_contig_idx];
  
  Assert( ins.IsFirstAnchor( first_anchor_contig.Id() ) );

  insertion_anchor first_anchor = ins.FirstAnchor();
  contig_range first_anchor_range = first_anchor.AnchorRange();

  int last_anchor_entry_idx = -1;
  for ( unsigned int ii = first_anchor_entry_idx; ii < entries_.size(); ++ii )
  {
    if ( entries_[ii].type == contig_type )
    {
      int contig_idx = entries_[ii].index;
      agp_contig& ctg = contigs_[ contig_idx ];

      if ( ins.IsLastAnchor( ctg.Id() ) )
      {
	last_anchor_entry_idx = ii;
	break;
      }
    }
  }

  if ( last_anchor_entry_idx == -1 )
  {
    cout << "Found only one anchor:" << endl; 
    cout << first_anchor << endl;
    cout << "for this insertion:" << endl;
    cout << "INSERT " << ins << endl;
    Assert( last_anchor_entry_idx > 0 );
  }

  agp_contig& last_anchor_contig = contigs_[entries_[last_anchor_entry_idx].index];
  insertion_anchor last_anchor = ins.LastAnchor();
  contig_range last_anchor_range = last_anchor.AnchorRange();

  // if the anchors are oriented the same way on the chromosome
  if ( first_anchor_contig.IsRC() == last_anchor_contig.IsRC() )
  {
    // but they are oriented differently with respect to the inserted contig
    if ( first_anchor.SameOrientation() != last_anchor.SameOrientation() )
    {
      // there's a problem
      cerr << "Mismatch of orientations for insertion:" << endl;
      cerr << ins << endl;
      Assert( false );
    }
  }

  // or if the anchors are not oriented the same way on the chromosome
  else
  {
    // but they are oriented the same way with respect to the inserted contig
    if ( first_anchor.SameOrientation() == last_anchor.SameOrientation() )
    {
      // there's a problem
      cerr << "Mismatch of orientations for insertion:" << endl;
      cerr << ins << endl;
      Assert( false );
    }
  }


  // if the anchor is reverse on the chromosome, trim its front
  if ( first_anchor.AnchorIsReversed() )
    first_anchor_contig.SetStart( first_anchor.AnchorRange().EndPos() + 1 );
  
  // if the anchor is forward on the chromosome, trim its back
  else 
    first_anchor_contig.SetStop( first_anchor.AnchorRange().BeginPos() - 1);
  
  vec<entry_pair>::iterator replacement_start_iter =
    entries_.begin() + first_anchor_entry_idx;

  // if the anchor has not been trimmed to nothing
  if ( first_anchor_contig.Extent() > 0 )
  {
    // keep it
    ++replacement_start_iter;
    
    // trim the new contig to the end of the match with the anchor
    if ( new_contig.IsRC() )
      new_contig.SetStop( first_anchor.ReplacementRange().EndPos() );
    else
      new_contig.SetStart( first_anchor.ReplacementRange().BeginPos() );
  }

  // else the anchor has been trimmed to nothing
  else
  {
    int prev_entry_idx = first_anchor_entry_idx - 1;
    if ( prev_entry_idx > 0 )
    {
      int replacement_overhang;
      contig_range replacement_range = first_anchor.ReplacementRange();
      if ( new_contig.IsRC() )
	replacement_overhang = 
	   new_contig.Length() - replacement_range.EndPos() - 1;
      else
	replacement_overhang = replacement_range.BeginPos();

      if ( entries_[prev_entry_idx].type == gap_type )
      {
	agp_gap& prev_gap = gaps_[ entries_[prev_entry_idx].index ];
	prev_gap.SetLength( max<int>( prev_gap.Length() - replacement_overhang,
				      minimum_gap_size ) );
      }
      else
      {
	int prev_ctg_trim_amount = -first_anchor_contig.Extent() + replacement_overhang;
	
	agp_contig& prev_ctg = contigs_[ entries_[prev_entry_idx].index ];

	cout << "While inserting " << new_contig.Name()
	     << " needed to trim " << prev_ctg.Name()
	     << " by " << prev_ctg_trim_amount
	     << " (" << -first_anchor_contig.Extent() << "+" << replacement_overhang << ")";

	if ( prev_ctg.Extent() <= prev_ctg_trim_amount )
	{
	  cout << ", but " << prev_ctg.Name() 
	       << " has only " << prev_ctg.Length() << " bases." << endl;
	  AssertGt( prev_ctg.Extent(), prev_ctg_trim_amount );
	}
	else
	  cout << "." << endl;

	if ( prev_ctg.IsRC() )
	  prev_ctg.SetStart( prev_ctg.Start() + prev_ctg_trim_amount );
	else
	  prev_ctg.SetStop( prev_ctg.Stop() - prev_ctg_trim_amount );
      }
    }
  }
	
  // if the anchor is forward on the chromsome, trim the front
  if ( ! last_anchor.AnchorIsReversed() )
    last_anchor_contig.SetStart( last_anchor.AnchorRange().EndPos() + 1 );
      
  // if the anchor is reverse on the chromosome, trim the back
  else
    last_anchor_contig.SetStop( last_anchor.AnchorRange().BeginPos() - 1);

  vec<entry_pair>::iterator replacement_stop_iter = 
    entries_.begin() + last_anchor_entry_idx + 1;

  // if the anchor has not been trimmed to nothing
  if ( last_anchor_contig.Extent() > 0 )
  {
    --replacement_stop_iter;
    
    if ( new_contig.IsRC() )
      new_contig.SetStart( last_anchor.ReplacementRange().BeginPos() );
    else
      new_contig.SetStop( last_anchor.ReplacementRange().EndPos() );
    
  }

  // else the anchor has been trimmed to nothing
  else
  {
    unsigned int next_entry_idx = last_anchor_entry_idx + 1;
    if ( next_entry_idx < entries_.size() )
    {
      int replacement_overhang;
      contig_range replacement_range = last_anchor.ReplacementRange();
      if ( new_contig.IsRC() )
	replacement_overhang = replacement_range.BeginPos();
      else
	replacement_overhang = 
	  new_contig.Length() - replacement_range.EndPos() - 1;

      if ( entries_[next_entry_idx].type == gap_type )
      {
	agp_gap& next_gap = gaps_[ entries_[next_entry_idx].index ];
	next_gap.SetLength( max<int>( next_gap.Length() - replacement_overhang,
				      minimum_gap_size ) );
      }
      else
      {
	int next_ctg_trim_amount = -last_anchor_contig.Extent() + replacement_overhang;
	
	agp_contig& next_ctg = contigs_[ entries_[next_entry_idx].index ];

	cout << "While inserting " << new_contig.Name()
	     << " needed to trim " << next_ctg.Name()
	     << " by " << next_ctg_trim_amount
	     << " (" << -last_anchor_contig.Extent() << "+" << replacement_overhang << ")";

	if ( next_ctg.Extent() <= next_ctg_trim_amount )
	{
	  cout << ", but " << next_ctg.Name() 
	       << " has only " << next_ctg.Length() << " bases." << endl;
	  AssertGt( next_ctg.Extent(), next_ctg_trim_amount );
	}
	else
	  cout << "." << endl;

	if ( next_ctg.IsRC() )
	  next_ctg.SetStop( next_ctg.Stop() - next_ctg_trim_amount );
	else
	  next_ctg.SetStart( next_ctg.Start() + next_ctg_trim_amount );
      }
    }
  }

  cout << new_contig.Name() << " (" 
       << new_contig.Start() << "-" << new_contig.Stop()
       << " of " << new_contig.Length() << ")"
       << " will replace " << distance( replacement_start_iter,
					replacement_stop_iter )
       << " entries (including gaps), including:" << endl;

  for ( vec<entry_pair>::iterator iter = replacement_start_iter;
	iter != replacement_stop_iter;
	++iter )
    if ( iter->type == contig_type )
      cout << contigs_[ iter->index ].Name() << endl;
  
  entries_.erase( replacement_start_iter, replacement_stop_iter );
  entries_.insert( replacement_start_iter, 
		   entry_pair( contig_type, contigs_.size() ) );

  contigs_.push_back( new_contig );
  
  return true;
}



void agp_chromosome::Load( const String &file_name )
{
  ForceAssert( IsRegularFile( file_name ) );

  ifstream in( file_name.c_str( ) );
  
  String a_line;
  vec<String> tokens;

  while ( in ) {
    getline( in, a_line );
    if ( !in )
      break;

    Tokenize( a_line, tokens );

    const String &last_t = tokens[ tokens.size( ) -1 ];
    
    // a_line is an agp_contig.
    if ( "+" == last_t || "-" == last_t ) {
      if ( "" == name_ )
	name_ = tokens[0].After( "chr" );

      agp_contig::contig_type type;

      type = TypeFromChar( tokens[4][0] );
      
      bool is_RC = ( "-" == tokens[8] ) ? true : false;
      
      agp_contig new_contig;
      
      new_contig.SetName( tokens[5] );
      new_contig.SetLength( -1 );
      new_contig.SetStart( tokens[6].Int( ) - 1 );
      new_contig.SetStop( tokens[7].Int( ) - 1 );
      new_contig.SetType( type );
      new_contig.SetRC( is_RC );

      this->AddContig( new_contig );
    }
    
    // a_line is an agp_gap.
    else if ( "yes" == last_t || "no" == last_t ) {
      if ( "" == name_ )
	name_ = tokens[0].After( "chr" );
      
      bool is_bridged = ( "yes" == tokens[7] ) ? true : false;

      int len = tokens[2].Int( ) - tokens[1].Int( ) + 1;

      agp_gap new_gap;
      
      new_gap.SetType( tokens[6] );
      new_gap.SetLength( len );
      new_gap.SetBridged( is_bridged );

      this->AddGap( new_gap );
    }
    
    // An error.
    else
      ForceAssert( 1 == 0 );
  }

}



void agp_chromosome::Print( ostream& out, const String *base_name )
{
  longlong base_position = 1;
  longlong unit_position = 1;

  for ( unsigned int ii = 0; ii < entries_.size(); ++ii )
  {
    agp_entry* entry;
    if ( entries_[ii].type == contig_type )
      entry = &(contigs_[entries_[ii].index]);
    else
      entry = &(gaps_[entries_[ii].index]);
    entry->Print( out, this->Name(), base_position, unit_position, base_name );
  }
}

