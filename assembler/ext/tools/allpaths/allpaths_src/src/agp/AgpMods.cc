// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#include <strstream>

#include "agp/AgpMods.h"

#include "system/System.h"

bool 
agp_mods::SuperShouldBeRemoved( int super_id ) const
{
  return ( binary_search( supers_to_remove_.begin(),
			  supers_to_remove_.end(),
			  super_id ) );
}

bool 
agp_mods::ContigShouldBeUnplaced( int contig_id ) const
{
  return ( binary_search( contigs_to_unplace_.begin(),
			  contigs_to_unplace_.end(),
			  contig_id ) );
}



void
agp_mods::Read( const String& filename )
{
  Ifstream( mods_stream, filename );

  String line_type;
  
  while ( mods_stream )
  {
    mods_stream >> line_type;

    if ( ! mods_stream )
      break;

    if ( line_type.Contains( "//", 0 ) )
    {
      String rest_of_line;
      getline( mods_stream, rest_of_line );
      continue;
    }

    if ( line_type == "INSERT" )
    {
      insertion new_insertion;
      mods_stream >> new_insertion;

      insertions_.push_back( new_insertion );
    }

    else if ( line_type == "DELETE" )
    {
      String deletion_type, line;
      getline( mods_stream, line );
      istrstream iline( line.c_str( ) );
      int super_to_delete;
      iline >> deletion_type;
      ForceAssertEq( deletion_type, "SUPERS" );
      while(1)
      {    iline >> super_to_delete;
           supers_to_remove_.push_back( super_to_delete );
           if ( !iline ) break;    }
    }

    else if ( line_type == "UNPLACE" )
    {
      String unplacement_type;
      int contig_to_unplace;
      mods_stream >> unplacement_type >> contig_to_unplace;

      AssertEq( unplacement_type, "CONTIG" );
      
      contigs_to_unplace_.push_back( contig_to_unplace );
    }

    else if ( line_type == "TRIM" )
    {
      contig_trim the_trim;
      mods_stream >> the_trim;

      trims_.push_back( the_trim );
    }      
    
    else
    {
      InputErr( "Unknown command " << line_type << " in " << filename << "." );
    }
  }

  sort( supers_to_remove_.begin(), supers_to_remove_.end() );
  sort( contigs_to_unplace_.begin(), contigs_to_unplace_.end() );
}
 
void 
agp_mods::Write( const String& filename ) const
{
  Ofstream( mods_stream, filename );

  for ( unsigned int i = 0; i < insertions_.size(); ++i )
    mods_stream << "INSERT " << insertions_[i] << "\n";

  for ( unsigned int i = 0; i < supers_to_remove_.size(); ++i )
    mods_stream << "DELETE SUPER " << supers_to_remove_[i] << "\n";

  for ( unsigned int i = 0; i < contigs_to_unplace_.size(); ++i )
    mods_stream << "UNPLACE CONTIG " << contigs_to_unplace_[i] << "\n";

  for ( unsigned int i = 0; i < trims_.size(); ++i )
    mods_stream << "TRIM " << trims_[i] << "\n";
}

ostream& 
operator<< ( ostream& out, const contig_range& range )
{
  out << range.begin_pos_ << "-" << range.end_pos_;

  return out;
}

istream& 
operator>> ( istream& in, insertion_anchor& anchor )
{
  String line, word;
  getline( in, line );

  // parse line of form:
  // 123-1435 contig 12341 @ 0-1312

  contig_range range;

  if ( !line.Contains( "-" ) ) cout << "confused by " << line << endl;

  ForceAssert( line.Contains( "-" ) );
  word = line.Before( "-" );
  line = line.After( "-" );

  Assert( word.IsInt() );
  range.SetBeginPos( word.Int() );
  
  Assert( line.Contains( " contig " ) );
  word = line.Before ( " contig " );
  line = line.After( " contig " );

  Assert( word.IsInt() );
  range.SetEndPos( word.Int() );

  anchor.SetReplacementRange( range );

  Assert( line.Contains( "@" ) );
  word = line.Before( "@" );
  line = line.After( "@" );

  Assert( word.IsInt() );
  anchor.anchor_id_ = word.Int();

  Assert( line.Contains( "-" ) );
  word = line.Before( "-" );
  line = line.After( "-" );

  Assert( word.IsInt() );
  range.SetBeginPos( word.Int() );
  
  word = line;

  Assert( word.IsInt() );
  range.SetEndPos( word.Int() );

  anchor.SetAnchorRange( range );

  return in;
}
	  
ostream& 
operator<< ( ostream& out, const insertion_anchor& anchor )
{
  out << anchor.replacement_range_
      << " contig " << anchor.anchor_id_ << "@" << anchor.anchor_range_ << endl;

  return out;
}

istream& 
operator>> ( istream& in, insertion& the_insertion )
{
  in >> the_insertion.id_;

  in >> the_insertion.first_anchor_;
  in >> the_insertion.last_anchor_;
  
  return in;
}

ostream&
operator<< ( ostream& out, const insertion& the_insertion )
{
  out << the_insertion.id_ << "\n";

  out << the_insertion.first_anchor_;
  out << the_insertion.last_anchor_;
  
  return out;
}

istream& 
operator>> ( istream& in, contig_trim& the_trim )
{
  String word_from;
  String end;
  
  in >> the_trim.amount_ >> word_from >> end;
  
  AssertEq( word_from, "FROM" );
  if ( end == "FRONT" )
  {
    the_trim.front_ = true;
  }
  else if ( end == "BACK" )
  {
    the_trim.front_ = false;
  }
  else
    InputErr( "Malformed trim command." );
  
  String word_of;
  String word_contig;
  
  in >> word_of >> word_contig >> the_trim.id_;
  
  AssertEq( word_of, "OF" );
  AssertEq( word_contig, "CONTIG" );

  return in;
}

ostream&
operator<< ( ostream& out, const contig_trim& the_trim )
{
  out << the_trim.Amount() << " FROM "
      << ( the_trim.TrimFront() ? "FRONT " : "BACK " )
      << "OF CONTIG " << the_trim.Id() << endl;

  return out;
}
