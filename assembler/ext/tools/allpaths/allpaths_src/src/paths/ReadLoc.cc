///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "CoreTools.h"
#include "PackAlign.h"
#include "PairsManager.h"
#include "PrintAlignment.h"
#include "Qualvector.h"
#include "VecUtilities.h"
#include "math/Functions.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "paths/ReadLoc.h"
#include "polymorphism/DumbCall.h"

vec< vec<int> > read_loc::lib_sep_, read_loc::lib_dev_;

void read_loc::LoadSepDev( const String& run_dir )
{    lib_sep_.resize(3), lib_dev_.resize(3);
     longlong nreads;
     vec<String> lib_names;
     if ( IsRegularFile( run_dir + "/frag_reads_filt_cpd.pairs" ) )
     {    ReadPairsManagerLibInfo( run_dir + "/frag_reads_filt_cpd.pairs",
               nreads, lib_names, lib_sep_[0], lib_dev_[0] );    }
     if ( IsRegularFile( run_dir + "/jump_reads_filt_cpd.pairs" ) )
     {    ReadPairsManagerLibInfo( run_dir + "/jump_reads_filt_cpd.pairs",
               nreads, lib_names, lib_sep_[1], lib_dev_[1] );    }
     else if ( IsRegularFile( run_dir + "/jump_reads_filt.pairs" ) )
     {    ReadPairsManagerLibInfo( run_dir + "/jump_reads_filt.pairs",
               nreads, lib_names, lib_sep_[1], lib_dev_[1] );    }
     if ( IsRegularFile( run_dir + "/long_jump_reads_filt.pairs" ) )
     {    ReadPairsManagerLibInfo( run_dir + "/long_jump_reads_filt.pairs",
               nreads, lib_names, lib_sep_[2], lib_dev_[2] );    }    }

void LoadLibNames( const String& run_dir, vec< vec<String> >& lib_names ) {
  lib_names.resize(3);
  longlong nreads;
  vec<int> lib_sep, lib_dev;
  if ( IsRegularFile( run_dir + "/frag_reads_filt_cpd.pairs" ) ) {
    ReadPairsManagerLibInfo( run_dir + "/frag_reads_filt_cpd.pairs",
			     nreads, lib_names[0], lib_sep, lib_dev );
  }
  if ( IsRegularFile( run_dir + "/jump_reads_filt_cpd.pairs" ) ) {
    ReadPairsManagerLibInfo( run_dir + "/jump_reads_filt_cpd.pairs",
			     nreads, lib_names[1], lib_sep, lib_dev );
  } else if ( IsRegularFile( run_dir + "/jump_reads_filt.pairs" ) ) {
    ReadPairsManagerLibInfo( run_dir + "/jump_reads_filt.pairs",
			     nreads, lib_names[1], lib_sep, lib_dev );
  }
  if ( IsRegularFile( run_dir + "/long_jump_reads_filt.pairs" ) ) {
    ReadPairsManagerLibInfo( run_dir + "/long_jump_reads_filt.pairs",
			     nreads, lib_names[2], lib_sep, lib_dev );
  }
}

int read_loc::Sep( ) const
{    int c = ReadClass( );
     if ( lib_sep_.size( ) == 0 || lib_sep_[c].size( ) == 0 )
     {    cout << "Sep called for " << ReadClassName( ) << ", but lib info "
               << "not loaded." << endl << "Abort." << endl;
          TracebackThisProcess( );    }
     return lib_sep_[c][ LibraryId( ) ];    }

int read_loc::Dev( ) const
{    int c = ReadClass( );
     if ( lib_dev_.size( ) == 0 || lib_dev_[c].size( ) == 0 )
     {    cout << "Dev called for " << ReadClassName( ) << ", but lib info "
               << "not loaded." << endl << "Abort." << endl;
          TracebackThisProcess( );    }
     return lib_dev_[c][ LibraryId( ) ];    }

String read_locs_packed_header = "This file represents a vector of read_locs\n";

void WriteReadLocs( const String& HEAD, const int32_t ntigs, 
     const vec<read_loc>& locs )
{    int fd = OpenForWrite( HEAD + ".readlocs" );
     WriteBytes( fd, &read_locs_packed_header[0], read_locs_packed_header.size( ) );
     size_t locs_size = locs.size( );
     WriteBytes( fd, &locs_size, sizeof(locs_size) );
     size_t bucket = 100000;
     vec<PackedStuff> locs_packed(bucket);
     for ( size_t i = 0; i < locs_size; i += bucket )
     {    size_t j;
          for ( j = 0; j < Min( bucket, locs_size - i ); j++ )
               locs_packed[j] = locs[i+j].PackWithoutContigId( );
          WriteBytes( fd, &locs_packed[0], sizeof(PackedStuff) * j );    }
     vec<int64_t> index( ntigs + 1 );
     uint64_t pos = 0;
     for ( int64_t u = 0; u <= (int64_t) ntigs; u++ )
     {    while( pos < locs.size( ) && locs[pos].ContigId( ) < u ) ++pos;
          index[u] = pos;    }
     WriteBytes( fd, &ntigs, sizeof(ntigs) );
     WriteBytes( fd, &index[0], (ntigs+1) * sizeof(int64_t) );
     Close(fd);    }

read_locs_on_disk::read_locs_on_disk( const String& HEAD, const String& run_dir )
: fr_((HEAD + ".readlocs").c_str())
{    read_loc dummy;
     dummy.LoadSepDev(run_dir);
     fr_.seek( read_locs_packed_header.size() );
     fr_.read( &locs_size_, sizeof(locs_size_) );
     fr_.seekRel( sizeof(PackedStuff) * locs_size_ );
     int32_t ntigs;
     fr_.read( &ntigs, sizeof(ntigs) );
     index_.resize(ntigs+1);
     fr_.read( &index_[0], (ntigs+1) * sizeof(int64_t) );    }

void read_locs_on_disk::LoadContig( const int tig, vec<read_loc>& locs )
{    ForceAssertGe( tig, 0 );
     ForceAssertLt( tig, index_.isize( ) );
     locs.resize( index_[tig+1] - index_[tig] );
     if ( locs.empty( ) ) return;
     fr_.seek( read_locs_packed_header.size( ) + sizeof(locs_size_)
          + sizeof(PackedStuff) * index_[tig] );
     vec<PackedStuff> locs_packed( locs.size( ) );
     fr_.read( &locs_packed[0], locs.size( ) * sizeof(PackedStuff) );
     for ( size_t i = 0; i < locs.size( ); i++ )
          locs[i] = read_loc( locs_packed[i], tig );    }

// Rewrite to the readlocs file. assuming the size of the locs for this contig12
// is not changed.
void read_locs_on_disk::WriteContig(const String& file,  const int tig, const vec<read_loc>& locs )
{    ForceAssertGe( tig, 0 );
     ForceAssertLt( tig, index_.isize( ) );
     ForceAssertEq( locs.size(),  size_t(index_[tig+1] - index_[tig]) );
     int fd = OpenForWrite(file );
     lseek( fd, read_locs_packed_header.size( ) + sizeof(locs_size_)
          + sizeof(PackedStuff) * index_[tig], SEEK_SET );
     vec<PackedStuff> locs_packed( locs.size( ) );
     for ( size_t i = 0; i < locs.size(); i++)
	  locs_packed[i] = locs[i].PackWithoutContigId( );
     WriteBytes( fd, &locs_packed[0], sizeof(PackedStuff) * locs.size() );
     Close(fd);    }

void read_locs_on_disk::LoadOne( const int64_t id, read_loc& rl )
{    ForceAssertGe( id, 0 );
     fr_.seek( read_locs_packed_header.size( ) + sizeof(locs_size_)
          + sizeof(PackedStuff) * id );
     PackedStuff rl_packed;
     fr_.read( &rl_packed, sizeof(PackedStuff) );
     int tig = upper_bound( index_.begin( ), index_.end( ), id ) 
          - index_.begin( ) - 1;
     rl = read_loc( rl_packed, tig );    }

String read_loc::ReadClassName( ) const
{    if ( ReadClass( ) == 0 ) return "frag";
     else if ( ReadClass( ) == 1 ) return "jump";
     else return "long";    }

ostream& operator<<( ostream& out, const read_loc& rl ) // OLD!
{    return out << rl.ReadClassName( ) << "." << rl.ReadId( ) 
          << " (len=" << rl.ReadLength( ) << ") --> " 
          << ( rl.Fw( ) ? "" : "rc " ) << "tig " << rl.ContigId( ) << " @ " 
          << rl.Start( ) << "\n";    }

void read_loc::PrintVisualLoc( const Bool visual_abbr, ostream& out,
     const String& read_name, const basevector& b, const basevector& t,
     const qualvector* q, const Bool partner ) const
{    align a;
     GetAlign( a, b, t, partner );
     out << read_name << ": " << a.pos1( ) << " - " << a.Pos1( ) << ", contig: " 
          << a.pos2( ) << " - " << a.Pos2( ) << " of " << t.size( ) << "\n ";
     ostringstream xout;
     if ( q == 0 ) PrintVisualAlignment( visual_abbr, xout, b, t, a );
     else PrintVisualAlignment( visual_abbr, xout, b, t, a, *q );
     String vis = xout.str( );
     for ( int i = 0; i < 80; i++ )
          vis.GlobalReplaceBy( " \n", "\n" );
     vis.GlobalReplaceBy( "\n\n\n", "\n\n" );
     vis.GlobalReplaceBy( "\n\n\n", "\n\n" );
     if ( vis.Contains( "\n\n", -1 ) ) vis.resize( vis.isize( ) - 1 );
     if ( vis.size( ) > 0 && vis[0] == '\n' ) vis = vis.After( "\n" );
     if ( vis.size( ) == 0 || vis[0] != '\n' ) out << "\n";
     vis.GlobalReplaceBy( "\n\n", "\n \n" );
     out << vis;    }

void PrintWithPairInfo( ostream& out, const vec<read_loc>& rl,
     const vec<superb>& scaffolds,
     const vec<int>& to_super, const vec<int>& to_super_pos,
     const vec<int>& to_super_posr, const String& run_dir, const String& sub_dir, 
     const String& ASSEMBLY, const Bool SHOW_ALIGNS, const Bool SHOW_PARTNER_ALIGNS,
     const Bool SHOW_QUALS_IN_ALIGNS )
{    int n = rl.size( );
     vec<int64_t> frag_read_ids, jump_read_ids, long_jump_read_ids;
     vec<int32_t> tig_ids;
     // Note that contig ids are getting truncated here by being stuffed into
     // int32_t instead of uint32_t.  The read ids are also getting truncated
     // but this can't possibly matter.
     vecbasevector frag_reads, jump_reads, long_jump_reads, tigs;
     vecqualvector frag_quals, jump_quals, long_jump_quals;
     read_locs_on_disk locs_file( sub_dir + "/" + ASSEMBLY, run_dir );
     if (SHOW_ALIGNS)
     {    for ( int i = 0; i < n; i++ )
          {    tig_ids.push_back( rl[i].ContigId( ) );
               if ( SHOW_PARTNER_ALIGNS && rl[i].PartnerPlaced( ) )
                    tig_ids.push_back( rl[i].PartnerContigId( ) );
               if ( rl[i].ReadClass( ) == 0 )
               {    frag_read_ids.push_back( rl[i].ReadId( ) );
                    if (SHOW_PARTNER_ALIGNS)
                         frag_read_ids.push_back( rl[i].PartnerReadId( ) );    }
               if ( rl[i].ReadClass( ) == 1 )
               {    jump_read_ids.push_back( rl[i].ReadId( ) );    
                    if (SHOW_PARTNER_ALIGNS)
                         jump_read_ids.push_back( rl[i].PartnerReadId( ) );    }
               if ( rl[i].ReadClass( ) == 2 )
               {    long_jump_read_ids.push_back( rl[i].ReadId( ) );
                    if (SHOW_PARTNER_ALIGNS)
                         long_jump_read_ids.push_back( 
                              rl[i].PartnerReadId( ) );    }    }
          UniqueSort(tig_ids);
          tigs.Read( sub_dir + "/" + ASSEMBLY + ".contigs.fastb", tig_ids );
          UniqueSort(frag_read_ids);
          UniqueSort(jump_read_ids);
          UniqueSort(long_jump_read_ids);
          if ( frag_read_ids.nonempty( ) )
               frag_reads.Read( run_dir + "/frag_reads_filt_cpd.fastb", frag_read_ids );
          if ( jump_read_ids.nonempty( ) )
               jump_reads.Read( run_dir + "/jump_reads_filt_cpd.fastb", jump_read_ids ); 
          if ( long_jump_read_ids.nonempty( ) )
          {    long_jump_reads.Read( run_dir + "/long_jump_reads_filt.fastb", 
                    long_jump_read_ids );    }
          if (SHOW_QUALS_IN_ALIGNS)
          {    if ( frag_read_ids.nonempty( ) )
               {    frag_quals.Read( run_dir + "/frag_reads_filt_cpd.qualb", 
                         frag_read_ids );    }
               if ( jump_read_ids.nonempty( ) )
               {    jump_quals.Read( run_dir + "/jump_reads_filt_cpd.qualb", 
                         jump_read_ids );    }
               if ( long_jump_read_ids.nonempty( ) )
               {    long_jump_quals.Read( run_dir + "/long_jump_reads_filt.qualb", 
                         long_jump_read_ids );    }    }    }
     double total_off = 0.0;
     int count_off = 0;
     for ( int i = 0; i < n; i++ )
     {    if (SHOW_ALIGNS)
          {    out << "=========================================================="
                    << "======================\n \n";    }
          uint32_t m1 = rl[i].ContigId( );
          int s1 = to_super[m1], p1 = to_super_pos[m1];
          out << rl[i].ReadClassName( ) << "[" << int(rl[i].LibraryId( )) << "]." 
               << rl[i].ReadId( ) << " (len=" << rl[i].ReadLength( ) << ") --> " 
               << ( rl[i].Fw( ) ? "fw " : "rc " ) << "tig " << m1 << " (s" 
               << to_super[m1] << "." << to_super_pos[m1] + 1 << "/" 
               << to_super_posr[m1] + 1 << ", len=" << scaffolds[s1].Len(p1)
               << ")" << " @ " << rl[i].Start( ) << "\n";
          if ( !rl[i].PartnerPlaced( ) ) 
          {    out << "(paired - " << rl[i].PartnerReadId( ) << ", unplaced)\n";
                    }
          else if ( rl[i].PartnerContigId( ) != m1 )
          {    int m2 = rl[i].PartnerContigId( );
               out << "(paired - " << rl[i].PartnerReadId( ) << " "
                    << ( rl[i].PartnerFw( ) ? "fw" : "rc" ) << " on contig " << m2
                    << " (s" << to_super[m2] << "." << to_super_pos[m2] + 1 
                    << "/" << to_super_posr[m2] + 1 << ")"
                    << " at position " << rl[i].PartnerStart( ) << "\n";    }
          else
          {    int delta = rl[i].SepDelta( );
               Bool catywampus = rl[i].Catywampus( );
               Bool backwards = rl[i].Backwards( );
               Bool degenerate = rl[i].Degenerate( );
               Bool first = ( rl[i].Start( ) < rl[i].PartnerStart( )
                    || ( rl[i].Start( ) == rl[i].PartnerStart( )
                         && rl[i].Fw( ) < rl[i].PartnerFw( ) )
                    || ( rl[i].Start( ) == rl[i].PartnerStart( )
                         && rl[i].Fw( ) == rl[i].PartnerFw( )
                         && rl[i].ReadId( ) < rl[i].PartnerReadId( ) ) );
               out << "(paired - " << rl[i].PartnerReadId( ) << " "
                    << ( first ? "below" : "above" ) << " "
                    << ( rl[i].PartnerFw( ) ? "fw" : "rc" ) << " @ "
                    << rl[i].PartnerStart( ) << ", ";
               if (catywampus) out << "catywampus";
               else if (backwards) out << "backwards";
               else if (degenerate) out << "degenerate";
               else 
               {    double offby = double(delta) / double( rl[i].Dev( ) );
                    if ( rl[i].ReadClass( ) == 1 )
                    {    if (first)
                         {    total_off += offby;
                              count_off += 1;    }
                         else 
                         {    total_off -= offby;
                              count_off -= 1;    }    }
                    out << "sep off " << delta << " = "
                        << std::fixed << setprecision(2)
                        << offby << " devs";
                    if ( rl[i].ReadClass( ) == 1 )
                    {    out << ", rolling mean = " 
                              << total_off / sqrt(double(count_off));    }    }
               out << ")\n";    }
          if (SHOW_ALIGNS)
          {    int npasses = ( SHOW_PARTNER_ALIGNS ? 2 : 1 );
               for ( int pass = 1; pass <= npasses; pass++ )
               {    uint64_t id 
                         = ( pass == 1 ? rl[i].ReadId( ) : rl[i].PartnerReadId( ) );
                    basevector b;
                    qualvector q;
                    if ( rl[i].ReadClass( ) == 0 )
                    {    b = frag_reads[ BinPosition( frag_read_ids, id ) ];
                         if (SHOW_QUALS_IN_ALIGNS)
                              q = frag_quals[ BinPosition(frag_read_ids, id) ];    }
                    if ( rl[i].ReadClass( ) == 1 )
                    {    b = jump_reads[ BinPosition( jump_read_ids, id ) ];
                         if (SHOW_QUALS_IN_ALIGNS)
                              q = jump_quals[ BinPosition( jump_read_ids, id ) ];   }
                    if ( rl[i].ReadClass( ) == 2 )
                    {    int64_t r = BinPosition( long_jump_read_ids, id );
                         b = long_jump_reads[r];
                         if (SHOW_QUALS_IN_ALIGNS) q = long_jump_quals[r];    }
                    if ( ( pass == 1 && rl[i].Rc( ) ) 
                         || ( pass == 2 && rl[i].PartnerRc( ) ) )
                    {    b.ReverseComplement( );
                         q.ReverseMe( );    }
                    if ( pass == 2 && !rl[i].PartnerPlaced( ) )
                    {    out << " \n";
                         b.Print( out, "partner" );
                         continue;    }
                    int m = ( 
                         pass == 1 ? rl[i].ContigId( ) : rl[i].PartnerContigId( ) );
                    const basevector& t = tigs[ BinPosition(tig_ids, m) ];
                    out << " \n";
                    rl[i].PrintVisualLoc( True, out, 
                         ( pass == 1 ? "read" : "partner read" ), b, t,
                         ( SHOW_QUALS_IN_ALIGNS ? &q : 0 ), 
                         ( pass == 2 ) );    }    }
          out << "\n";    }    }

// TODO: Implement scaffold alignments
void PrintSAMAligns( ostream& sam, const vec<read_loc>& rl,
		     const vec<int>& to_super, const vec<int>& to_super_pos,
		     const vec< vec<int> >& to_super_offset,
		     const String& run_dir, const String& sub_dir,
		     const String& ASSEMBLY ) {
  // values for the SAM flag field
  int const IS_PAIRED = 0x0001;
  int const IS_PROPERLY_PAIRED = 0x0002;
  int const NOT_MAPPED = 0x0004;
  int const MATE_NOT_MAPPED = 0x0008;
  int const IS_REVERSED = 0x0010;
  int const MATE_IS_REVERSED = 0x0020;
  int const FIRST_READ = 0x0040;
  int const SECOND_READ = 0x0080;
  int const NOT_PRIMARY = 0x0100;
  int const NOT_PF = 0x0200;

  int n = rl.size( );
  if ( n < 1 ) return;   // early exit - no read_locs.

  vec<int64_t> frag_read_ids, jump_read_ids, long_jump_read_ids;
  vec<int32_t> tig_ids;
  // Note that contig ids are getting truncated here by being stuffed into
  // int32_t instead of uint32_t.  The read ids are also getting truncated
  // but this can't possibly matter.
  vecbasevector frag_reads, jump_reads, long_jump_reads, tigs;
  vecqualvector frag_quals, jump_quals, long_jump_quals;
  vec< vec<String> > lib_names;
  read_locs_on_disk locs_file( sub_dir + "/" + ASSEMBLY, run_dir );
  LoadLibNames( run_dir, lib_names );
  for ( int i = 0; i < n; i++ ) {
    tig_ids.push_back( rl[i].ContigId( ) );
    if (rl[i].PartnerPlaced( ) )
      tig_ids.push_back( rl[i].PartnerContigId( ) );
    if ( rl[i].ReadClass( ) == 0 ) {
      frag_read_ids.push_back( rl[i].ReadId( ) );
      frag_read_ids.push_back( rl[i].PartnerReadId( ) );
    }
    if ( rl[i].ReadClass( ) == 1 ) {
      jump_read_ids.push_back( rl[i].ReadId( ) );
      jump_read_ids.push_back( rl[i].PartnerReadId( ) );
    }
    if ( rl[i].ReadClass( ) == 2 ) {
      long_jump_read_ids.push_back( rl[i].ReadId( ) );
      long_jump_read_ids.push_back(rl[i].PartnerReadId( ) );
    }
  }
  UniqueSort(tig_ids);
  tigs.Read( sub_dir + "/" + ASSEMBLY + ".contigs.fastb", tig_ids );
  UniqueSort(frag_read_ids);
  UniqueSort(jump_read_ids);
  UniqueSort(long_jump_read_ids);
  if ( frag_read_ids.nonempty( ) ) {
    frag_reads.Read( run_dir + "/frag_reads_filt_cpd.fastb", frag_read_ids );
    frag_quals.Read( run_dir + "/frag_reads_filt_cpd.qualb", frag_read_ids );
  }
  if ( jump_read_ids.nonempty( ) ) {
    String jump_head;
    if ( IsRegularFile( run_dir + "/jump_reads_filt_cpd.fastb" ) ) {
      jump_head = "/jump_reads_filt_cpd.";
    } else {
      jump_head = "/jump_reads_filt.";
    }
    jump_reads.Read( run_dir + jump_head + "fastb", jump_read_ids );
    jump_quals.Read( run_dir + jump_head + "qualb", jump_read_ids );
  }
  if ( long_jump_read_ids.nonempty( ) ) {
    long_jump_reads.Read( run_dir + "/long_jump_reads_filt.fastb", 
			  long_jump_read_ids );
    long_jump_quals.Read( run_dir + "/long_jump_reads_filt.qualb",
			  long_jump_read_ids );
  }

  for ( int i = 0; i < n; i++ ) {
    // From the SAM 1.4 Specification:
    // For sequencing data, fragments are indexed by the order in which they are
    // sequenced. For fragments of an assembled sequence, they are indexed by
    // the order of the leftmost coordinate of the assembled sequence.
    if ( ( !rl[i].PartnerPlaced( ) )
	 || ( rl[i].ContigId( ) == rl[i].PartnerContigId( )
	      && ( rl[i].Start( ) < rl[i].PartnerStart( )
		   || ( rl[i].Start( ) == rl[i].PartnerStart( ) 
			&& rl[i].ReadId( ) < rl[i].PartnerReadId( ) ) ) ) 
	 || rl[i].ContigId( ) < rl[i].PartnerContigId( ) ) {
      String rg_name = rl[i].ReadClassName( ) + "." +
	lib_names[rl[i].ReadClass( )][rl[i].LibraryId( )];
      String sam_qname = rg_name + "." +
	ToString( (uint)rl[i].PseudoPairId( ) );
      int sam_flag = 0, mate_sam_flag = 0;
      String sam_rname = "contig_" + ToString( (uint)rl[i].ContigId( ) );
      String mate_sam_rname;
      int sam_pos = 0;
      int sam_mapq = 255, mate_sam_mapq = 255;
      String sam_cigar = "", mate_sam_cigar = "";
      String sam_rnext = "=", mate_sam_rnext = "=";
      int sam_pnext = 0;
      int sam_tlen = 0;
      String sam_seq, mate_sam_seq;
      String sam_qual, mate_sam_qual;
      String sam_opt_fields = "RG:Z:" + rg_name;
      uint32_t m1 = rl[i].ContigId( );
      int s1 = to_super[m1], p1 = to_super_pos[m1];
      int o1 = to_super_offset[s1][p1];
      sam_flag |= IS_PAIRED;
      sam_flag |= FIRST_READ;
      mate_sam_flag |= IS_PAIRED;
      mate_sam_flag |= SECOND_READ;
      align r_align;
      for ( int pass = 1; pass <= 2; pass++ ) {
	uint64_t id = ( pass == 1 ? rl[i].ReadId( ) : rl[i].PartnerReadId( ) );
	basevector b;
	qualvector q;
	if ( rl[i].ReadClass( ) == 0 ) {
	  b = frag_reads[ BinPosition( frag_read_ids, id ) ];
	  q = frag_quals[ BinPosition( frag_read_ids, id ) ];
	} else if ( rl[i].ReadClass( ) == 1 ) {
	  b = jump_reads[ BinPosition( jump_read_ids, id ) ];
	  q = jump_quals[ BinPosition( jump_read_ids, id ) ];
	} else if ( rl[i].ReadClass( ) == 2 ) {
	  b = long_jump_reads[ BinPosition( long_jump_read_ids, id ) ];
	  q = long_jump_quals[ BinPosition( long_jump_read_ids, id ) ];
	}
	if ( ( pass == 1 && rl[i].Rc( ) ) 
	     || ( pass == 2 && rl[i].PartnerRc( ) ) ) {
	  b.ReverseComplement( );
	  q.ReverseMe( );
	  if ( pass == 1 ) {
	    sam_flag |= IS_REVERSED;
	    mate_sam_flag |= MATE_IS_REVERSED;
	  } else {
	    sam_flag |= MATE_IS_REVERSED;
	    mate_sam_flag |= IS_REVERSED;
	  }
	}
	if ( pass == 1 ) {
	  sam_seq = b.ToString( );
	  sam_qual.resize( sam_seq.size( ) );
	  for ( qvec::size_type j = 0; j < q.size( ); j++ ) {
	    sam_qual[j] = q[j] + 33;
	  }
	} else {
	  mate_sam_seq = b.ToString( );
	  mate_sam_qual.resize( mate_sam_seq.size( ) );
	  for ( qvec::size_type j = 0; j < q.size( ); j++ ) {
	    mate_sam_qual[j] = q[j] + 33;
	  }
	}
	if ( pass == 1 || rl[i].PartnerPlaced( ) ) {
	  int m = ( pass == 1 ? rl[i].ContigId( ) : rl[i].PartnerContigId( ) );
	  const basevector& t = tigs[ BinPosition(tig_ids, m) ];
	  align a;
	  rl[i].GetAlign( a, b, t, ( pass == 2 ) );
	  if ( pass == 1 ) {
	    r_align = a;
	    int bases_aligned = 0;
	    int nblocks = a.Nblocks( );
	    sam_pos = a.StartOnTarget( ) + 1;
	    if ( a.StartOnQuery( ) > 0 ) {
	      sam_cigar += ToString( a.StartOnQuery( ) ) + "S";
	      bases_aligned += a.StartOnQuery( );
	    }
	    if ( nblocks > 0 ) {
	      sam_cigar += ToString( a.Lengths( 0 ) ) + "M";
	      bases_aligned += a.Lengths( 0 );
	      for (int j = 1; j < nblocks; j++ ) {
		int a_gap = a.Gaps( j );
		if ( a_gap < 0 ) {
		  sam_cigar += ToString( -a_gap ) + "I";
		  bases_aligned -= a_gap;
		} else {
		  sam_cigar += ToString( a_gap ) + "D";
		}
		sam_cigar += ToString( a.Lengths( j ) ) + "M";
		bases_aligned += a.Lengths( j );
	      }
	    }
	    if ( sam_seq.size( ) != (uint)bases_aligned ) {
	      ForceAssertLt( bases_aligned, (int)sam_seq.size( ) );
	      sam_cigar += ToString( sam_seq.size( ) - bases_aligned ) + "S";
	    }
	  } else {
	    int bases_aligned = 0;
	    int nblocks = a.Nblocks( );
	    sam_pnext = a.StartOnTarget( ) + 1;
	    if ( a.StartOnQuery( ) > 0 ) {
	      mate_sam_cigar += ToString( a.StartOnQuery( ) ) + "S";
	      bases_aligned += a.StartOnQuery( );
	    }
	    if ( nblocks > 0 ) {
	      mate_sam_cigar += ToString( a.Lengths( 0 ) ) + "M";
	      bases_aligned += a.Lengths( 0 );
	      for (int j = 1; j < nblocks; j++ ) {
		int a_gap = a.Gaps( j );
		if ( a_gap < 0 ) {
		  mate_sam_cigar += ToString( -a_gap ) + "I";
		  bases_aligned -= a_gap;
		} else {
		  mate_sam_cigar += ToString( a_gap ) + "D";
		}
		mate_sam_cigar += ToString( a.Lengths( j ) ) + "M";
		bases_aligned += a.Lengths( j );
	      }
	    }
	    if ( mate_sam_seq.size( ) != (uint)bases_aligned ) {
	      ForceAssertLt( bases_aligned, (int)mate_sam_seq.size( ) );
	      mate_sam_cigar += ToString( mate_sam_seq.size( ) - bases_aligned )
		+ "S";
	    }
	    if  ( rl[i].PartnerContigId( ) == m1 ) {
	      int r_end = r_align.EndOnTarget( ) + 1;
	      int m_end = a.EndOnTarget( ) + 1;
	      sam_tlen = Max( m_end, r_end ) - Min( sam_pnext, sam_pos ) + 1;
	      if ( sam_pos > sam_pnext ||
		   ( sam_pos == sam_pnext &&
		     ( r_end > m_end ||
		       ( r_end == m_end && rl[i].Rc( ) &&
			 rl[i].PartnerFw( ) ) ) ) ) {
		sam_flag ^= FIRST_READ;
		sam_flag |= SECOND_READ;
		mate_sam_flag ^= SECOND_READ;
		mate_sam_flag |= FIRST_READ;
		sam_tlen *= -1;
	      }
	      int delta = rl[i].SepDelta( );
	      double offby = double(delta) / double( rl[i].Dev( ) );
	      mate_sam_rname = sam_rname;

	      // Consider alignment "proper" if it's within 5 sd
	      if ( offby * offby <= 25.0 ) {
		sam_flag |= IS_PROPERLY_PAIRED;
		mate_sam_flag |= IS_PROPERLY_PAIRED;
	      }
	    } else {
	      mate_sam_rname = sam_rnext = "contig_"
		+ ToString( (uint)rl[i].PartnerContigId( ) );
	      mate_sam_rnext = sam_rname;
	    }
	  }
	} else {
	  sam_flag |= MATE_NOT_MAPPED;
	  mate_sam_flag |= NOT_MAPPED;
	  mate_sam_rname = sam_rname;
	  sam_pnext = sam_pos;
	  mate_sam_cigar = "*";
	  mate_sam_mapq = 0;
	}
      }
      sam << sam_qname << "\t" << sam_flag << "\t" << sam_rname << "\t"
	  << sam_pos << "\t" << sam_mapq << "\t" << sam_cigar << "\t"
	  << sam_rnext << "\t" << sam_pnext << "\t" << sam_tlen << "\t"
	  << sam_seq << "\t" << sam_qual << "\t" << sam_opt_fields << "\n"
	  << sam_qname << "\t" << mate_sam_flag << "\t" << mate_sam_rname
	  << "\t" << sam_pnext << "\t" << mate_sam_mapq << "\t"
	  << mate_sam_cigar << "\t" << mate_sam_rnext << "\t" << sam_pos
	  << "\t" << ( -1 * sam_tlen ) << "\t" << mate_sam_seq << "\t"
	  << mate_sam_qual << "\t" << sam_opt_fields << "\n";
    }
  }
}


    read_loc::read_loc( const align& a, const uint64_t read_id, 
         const uint32_t contig_id, const bool fw_on_contig, 
         const uint8_t read_class, const uint8_t library_id, 
         const uint16_t read_length ) :
    read_id_(read_id), contig_id_(contig_id), read_length_(read_length),
    library_id_(library_id), fw_on_contig_(fw_on_contig), read_class_(read_class)
  { AssertLt(read_id,1ul<<33);
    AssertLt(read_class,1ul<<4);
    AssertLt(library_id,1ul<<8);
    AssertLt(read_length,1ul<<16);
    int p1 = a.pos1( ), p2 = a.pos2( );
    int pos_low = p2 - p1, pos_high = p2 - p1;
    for ( int j = 0; j < a.Nblocks( ); j++ )
    {    if ( a.Gaps(j) > 0 ) 
         {    p2 += a.Gaps(j);
              pos_high = Max( pos_high, p2 - p1 );    }
         if ( a.Gaps(j) < 0 ) 
         {    p1 -= a.Gaps(j);
              pos_low = Min( pos_low, p2 - p1 );    }    }
    pos_on_contig_ = ( pos_low + pos_high ) / 2;
    bandwidth_ = ( pos_high - pos_low + 1 ) / 2;
    partner_placed_ = false;
    partner_read_id_ = 0;
    partner_contig_id_ = 0;
    partner_pos_on_contig_ = 0;
    partner_fw_on_contig_ = false;
    partner_read_length_ = 0;

    // The following compensates for an undiagnosed bug in SmithWatBandedA.

    if ( bandwidth_ > 0 ) bandwidth_++;

    AssertLt( pos_on_contig_, 1l<<34 );
  }

void AddToPileup( const read_loc& rl, const basevector& b, const basevector& tig,
     vec<dumbcall>& calls )
{    align a;
     rl.GetAlign( a, b, tig );
     int p1 = a.pos1( ), p2 = a.pos2( );
     for ( int j = 0; j < a.Nblocks( ); j++ )
     {    if ( a.Gaps(j) > 0 )
          {    for ( int u = 0; u < a.Gaps(j); u++ )
                    calls[p2+u].base[4]++;
               p2 += a.Gaps(j);    }
          if ( a.Gaps(j) < 0 )
          {    for ( int u = 0; u < -a.Gaps(j); u++ )
                    calls[p2].base[5]++;
               p1 -= a.Gaps(j);    }
          for ( int x = 0; x < a.Lengths(j); x++ )
          {    calls[p2].base[ b[p1] ]++;
               ++p1; ++p2;    }    }    }

void Pileup( const basevector& tig, const vec<read_loc>& rl,
     const vecbasevector& frag_reads, const vecbasevector& jump_reads,
     const vecbasevector& long_jump_reads, vec<dumbcall>& calls )
{    calls.resize_and_set( tig.size( ), dumbcall( ) );
     for ( size_t i = 0; i < rl.size( ); i++ )
     {    uint64_t id = rl[i].ReadId( );
          basevector b;
          if ( rl[i].ReadClass( ) == 0 ) b = frag_reads[id];
          if ( rl[i].ReadClass( ) == 1 ) b = jump_reads[id];
          if ( rl[i].ReadClass( ) == 2 ) b = long_jump_reads[id];
          if ( rl[i].Rc( ) ) b.ReverseComplement( );
          AddToPileup( rl[i], b, tig, calls );    }    }

void Pileup( const basevector& tig, const vec<read_loc>& rl,
	     const String& run_dir, vec<dumbcall>& calls, const Bool get_frags,
	     const Bool get_jumps, const Bool get_longs ){    
  calls.resize_and_set( tig.size( ), dumbcall( ) );
  int n = rl.size( );
  vec<int64_t> frag_read_ids, jump_read_ids, long_jump_read_ids; 
  vecbasevector frag_reads, jump_reads, long_jump_reads;
  for ( int i = 0; i < n; i++ ){    
    if ( rl[i].ReadClass( ) == 0 && get_frags ) 
      frag_read_ids.push_back( rl[i].ReadId( ) );
    else if ( rl[i].ReadClass( ) == 1 && get_jumps ) 
      jump_read_ids.push_back( rl[i].ReadId( ) );
    else if ( rl[i].ReadClass( ) == 2 && get_longs ) 
      long_jump_read_ids.push_back( rl[i].ReadId( ) );    
  }
  UniqueSort(frag_read_ids);
  UniqueSort(jump_read_ids);
  UniqueSort(long_jump_read_ids);
  if ( frag_read_ids.nonempty( ) )
    frag_reads.Read( run_dir + "/frag_reads_filt_cpd.fastb", frag_read_ids );
 
  if ( jump_read_ids.nonempty( ) )
    jump_reads.Read( run_dir + "/jump_reads_filt_cpd.fastb", jump_read_ids );

  if ( long_jump_read_ids.nonempty( ) )
    long_jump_reads.Read( run_dir + "/long_jump_reads_filt.fastb", 
			  long_jump_read_ids );
  for ( size_t i = 0; i < rl.size( ); i++ ){    
    int id = rl[i].ReadId( );
    basevector b;
    if ( rl[i].ReadClass( ) == 0 && get_frags )
      b = frag_reads[ BinPosition( frag_read_ids, id ) ];
    else if ( rl[i].ReadClass( ) == 1 && get_jumps )
      b = jump_reads[ BinPosition( jump_read_ids, id ) ];
    else if ( rl[i].ReadClass( ) == 2 && get_longs )
      b = long_jump_reads[ BinPosition( long_jump_read_ids, id ) ];

    if ( b.size() > 0 ){
      if ( rl[i].Rc( ) ) b.ReverseComplement( );
      AddToPileup( rl[i], b, tig, calls );
    }    
  }    
}

void Pileup( const basevector& tig, const vec<read_loc>& rl,
     const String& run_dir, vec<dumbcall>& calls ){    
  Pileup( tig, rl, run_dir, calls, True, True, True );
}

void read_loc::GetAlign( align& a, const basevector& b, const basevector& t,
     const Bool partner ) const
{    int errors;
     if ( !partner ) SmithWatBandedA( b, t, -Start( ), Bandwidth( ), a, errors );
     else
     {    SmithWatBandedA( b, t, -PartnerStart( ), 
               PartnerBandwidth( ), a, errors );    }    }
