///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Compute physical coverage of selected libraries.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "MainTools.h"
#include "PairsManager.h"
#include "Superb.h"
#include "VecUtilities.h"
#include "paths/ReadLoc.h"
#include "reporting/PerfStat.h"

class CLibPhysical {

public:

  CLibPhysical ( ) : n_pairs_ ( 0 ), tot_bases_ ( 0 ), tot_cov_ ( 0.0 ) { }
  
  int64_t NPairs( ) const { return n_pairs_; }
  int64_t TotBases( ) const { return tot_bases_; }
  double TotCov( ) const { return tot_cov_; }
  double PCov( ) const { return tot_bases_ == 0 ? 0.0 : tot_cov_ / double( tot_bases_ ); }
  
  void Add( int64_t n_pairs, int64_t tot_bases, double tot_cov ) {
    n_pairs_ += n_pairs;
    tot_bases_ += tot_bases;
    tot_cov_ += tot_cov;
  }
  
  
private:

  int64_t n_pairs_;
  int64_t tot_bases_;
  double tot_cov_;

};

int main(int argc, char *argv[])
{
     RunTime( );
  
     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String_OrDefault(SUBDIR, "test");
     CommandArgument_String_OrDefault(ASSEMBLY, "linear_scaffolds0.clean.patched");
     CommandArgument_String_OrDefault(S, "all");
     CommandArgument_String_OrDefault_Doc(READLOCS_PREFIX, "",
          "if specified, file extension is .READLOCS_PREFIX.readlocs "
          "instead of .readlocs");
     CommandArgument_String_OrDefault(LIB_TYPES, "{frag,jump,long_jump}");
     CommandArgument_Double_OrDefault_Doc( MAX_DEV, 5.0,
	  "only consider inserts with length within the given library separation "
	  "plus or minus MAX_DEV deviations" );
     CommandArgument_Double_OrDefault_Doc( EXCLUDE_SEP, 1.0,
	  "for a given library, only consider the physical coverage of the portion of scaffold "
	  "excluding the first and last EXCLUDE_SEP * lib_sep bases" );
     CommandArgument_Bool_OrDefault_Doc( ARCHIVE, False,
          "if true, save output as ../<ASSEMBLY>.LibCoverage.out" );
     CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
	   "Number of threads to use (use all available processors if set to 0)");
     EndCommandArguments;

     // Thread control

     NUM_THREADS = configNumThreads(NUM_THREADS);
     omp_set_num_threads( NUM_THREADS );

     // Define directories.

     String data_dir = PRE + "/" + DATA, run_dir = PRE + "/" + DATA + "/" + RUN;
     String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR; 

     // Output stream.
     
     ofstream archive_out;
     if ( ARCHIVE ) {
       String archive_file =  sub_dir + "/" + ASSEMBLY + ".LibCoverage.out";
       cout << "Sending output to " << archive_file << "\n" << endl;
       archive_out.open( archive_file.c_str( ) );
       PrintCommandPretty( archive_out );
     }
     ostream &out = ARCHIVE ? archive_out : * (ostream *) &cout;
     
     // Load scaffolds.

     int ntigs 
          = MastervecFileObjectCount( sub_dir + "/" + ASSEMBLY + ".contigs.fastb" );
     vec<superb> scaffolds;
     ReadSuperbs( sub_dir + "/" + ASSEMBLY + ".superb", scaffolds );
     vec<int> to_super( ntigs, -1 ), to_super_pos( ntigs, -1 );
     for ( int i = 0; i < scaffolds.isize( ); i++ )
     {    int n = scaffolds[i].Ntigs( );
          for ( int j = 0; j < n; j++ )
          {    to_super[ scaffolds[i].Tig(j) ] = i;
               to_super_pos[ scaffolds[i].Tig(j) ] = j;    }    }

     // Load pairing info.
     
     vec<Bool> lib_types( 3, False );
     vec<String> lib_types_s;
     ParseStringSet( LIB_TYPES, lib_types_s );
     for ( int i = 0; i < lib_types_s.isize( ); i++ )
     {    if ( lib_types_s[i] == "frag" ) lib_types[0] = True;
          else if ( lib_types_s[i] == "jump" ) lib_types[1] = True;
          else if ( lib_types_s[i] == "long_jump" ) lib_types[2] = True;
          else
          {    out << "Illegal LIB_TYPES." << endl;
               return 1;    }    }
     
     vec<String> pairs_files( 3, "" );
     for (int pass=0; pass<3; pass++) {
       String type_name = pass == 0 ? "frag" : ( pass == 1 ? "jump" : "long_jump" );
       type_name += "_reads_filt" + String( pass != 2 ? "_cpd" : "" ) + ".pairs";
       String file_name = run_dir + "/" + type_name;
       if ( ! IsRegularFile( file_name ) ) {
	 cout << "Unable to find pairs file: " << file_name << endl;
	 lib_types[pass] = False;   // Note! Reset lib_types if pairs file not found.
	 continue;
       }
       pairs_files[pass] = file_name;
     }

     vec<PairsManager> pairs;   // frag, jump, and long_jump
     for (int pass=0; pass<3; pass++) {
       if ( ! lib_types[pass] ) {
	 pairs.push_back( PairsManager( ) );
	 continue;
       }
       String type_name = pass == 0 ? "frag" : ( pass == 1 ? "jump" : "long_jump" );
       out << Date( ) << ": loading pairing info for " << type_name << endl;
       pairs.push_back( PairsManager( pairs_files[pass] ) );
     }

     out << Date( ) << ": found " << pairs[0].nLibraries() + pairs[1].nLibraries() 
       + pairs[2].nLibraries() << " libraries" << endl;
     
     // Valid windows of separations for various libraries (NB: separations).

     vec< vec< pair<int,int> > > valid_win( 3 );   // frag, jump, and long_jump
     for (int pass=0; pass<3; pass++) {
       if ( ! lib_types[pass] ) continue;
       valid_win[pass].reserve( pairs[pass].nLibraries( ) );
       for (size_t ii=0; ii<pairs[pass].nLibraries( ); ii++) {
	 int sep = pairs[pass].getLibrarySep( ii );
	 int dev = pairs[pass].getLibrarySD( ii );
	 int radius = int( MAX_DEV * double( dev ) );
	 valid_win[pass].push_back( make_pair( sep - radius, sep + radius ) );
       }
     }

     // Exclude lengths.

     vec< vec<int> > excludes( 3 );   // frag, jump, and long_jump
     for (int pass=0; pass<3; pass++) {
       if ( ! lib_types[pass] ) continue;
       excludes[pass].reserve( pairs[pass].nLibraries( ) );
       for (size_t ii=0; ii<pairs[pass].nLibraries( ); ii++) {
	 int sep = pairs[pass].getLibrarySep( ii );
	 int radius = int( EXCLUDE_SEP * double( sep ) );
	 excludes[pass].push_back( radius );
       }
     }
     
     // Raw lib data.

     vec< vec<CLibPhysical> > libdata( 3 );   // frag, jump, and long_jump
     vec< vec<size_t> > assembled_len( 3 );   // frag, jump, and long_jump
     vec< vec<size_t> > assembled_count( 3 );   // frag, jump, and long_jump
     for (int pass=0; pass<3; pass++) {
       if ( ! lib_types[pass] ) continue;
       int nlibs = pairs[pass].nLibraries( );
       libdata[pass].resize( nlibs );
       assembled_len[pass].resize( nlibs );
       assembled_count[pass].resize( nlibs );
     }

     out << Date( ) << ": computing physical coverage" << endl;

     // Compute physical coverage.

     vec<int> sids;
     if ( S != "all" ) ParseIntSet( S, sids );
     else
     {    for ( int j = 0; j < scaffolds.isize( ); j++ )
               sids.push_back(j);    }
     String head = sub_dir + "/" + ASSEMBLY;
     if ( READLOCS_PREFIX != "" ) head += "." + READLOCS_PREFIX;
     read_locs_on_disk locs_file( head, run_dir );
     #pragma omp parallel for
     for ( int z = 0; z < sids.isize( ); z++ )
     {    int s1 = sids[z];
          const superb& S1 = scaffolds[s1];
          vec< vec < vec< vec<int> > > > cov( 3 );
	  for (int pass=0; pass<3; pass++) {
	    if ( ! lib_types[pass] ) continue;
	    cov[pass].resize( pairs[pass].nLibraries( ) );
	    for (size_t jj=0; jj<pairs[pass].nLibraries( ); jj++) {
	      cov[pass][jj].resize(  S1.Ntigs( ) );
	      for ( int p1 = 0; p1 < S1.Ntigs( ); p1++ )
		cov[pass][jj][p1].resize( S1.Len(p1), 0 );	      
	    }
	  }
          double cov_frac 
               = double( S1.ReducedLength( ) ) / double( S1.FullLength( ) );
          vec< vec<int64_t> > nfosmids_this( 3 );
	  for (int pass=0; pass<3; pass++) {
	    if ( ! lib_types[pass] ) continue;
	    nfosmids_this[pass].resize( pairs[pass].nLibraries( ), 0 );
	  }
          for ( int p1 = 0; p1 < S1.Ntigs( ); p1++ )
          {    int m1 = S1.Tig(p1);
               vec<read_loc> locs;
               #pragma omp critical
               locs_file.LoadContig( m1, locs );
               for ( int j = 0; j < locs.isize( ); j++ )
               {    const read_loc& rl = locs[j];
		    int read_class = rl.ReadClass( );
		    if ( !lib_types[read_class] ) continue;
		    int lib_id = rl.LibraryId( );
		    assembled_len[read_class][lib_id] += size_t( rl.ReadLength( ) );
		    assembled_count[read_class][lib_id]++;
                    if ( !rl.PartnerPlaced( ) ) continue;
                    if ( !rl.Fw( ) || !rl.PartnerRc( ) ) continue;
                    int x1 = rl.Start( );
                    int y1 = S1.SubSuperLength( 0, p1 - 1 ) + x1;
                    int m2 = rl.PartnerContigId( );
                    int s2 = to_super[m2], p2 = to_super_pos[m2];
                    if ( s2 != s1 ) continue;
                    int x2 = rl.PartnerStart( );
                    const superb& S2 = scaffolds[s2];
                    int y2 = S2.SubSuperLength( 0, p2 - 1 ) + x2;
                    int len =  y2 + rl.PartnerReadLength( ) - y1;
		    int insert_sep = len - rl.ReadLength( ) - rl.PartnerReadLength( );
		    if ( insert_sep < valid_win[read_class][lib_id].first ) continue;
		    if ( insert_sep >= valid_win[read_class][lib_id].second ) continue;
                    nfosmids_this[read_class][lib_id]++;
                    if ( p1 == p2 )
                    {    for ( int u = x1; u < x2 + rl.PartnerReadLength( ); u++ )
                         {    if ( u >= 0 && u < cov[read_class][lib_id][p1].isize( ) ) 
                                   cov[read_class][lib_id][p1][u]++;     }    }
                    else
                    {    for ( int u = x1; u < cov[read_class][lib_id][p1].isize( ); u++ )
                              if ( u >= 0 ) cov[read_class][lib_id][p1][u]++;
                         for ( int p = p1 + 1; p < p2; p++ )
                         {    for ( int u = 0; u < cov[read_class][lib_id][p].isize( ); u++ )
                                   cov[read_class][lib_id][p][u]++;    }
                         for ( int u = 0; u < x2 + rl.PartnerReadLength( ); u++ )
                         {    if ( u < cov[read_class][lib_id][p2].isize( ) ) 
                                   cov[read_class][lib_id][p2][u]++;    }    }    }    }

          int N = S1.FullLength( );

	  for (int pass=0; pass<3; pass++) {
	    if ( ! lib_types[pass] ) continue;

	    for (size_t lib_id=0; lib_id<pairs[pass].nLibraries( ); lib_id++) {
	      int64_t pos = 0, total_bases_this = 0, total_cov_this = 0;
	      int exclude = excludes[pass][lib_id];
	      for ( int p = 0; p < S1.Ntigs( ); p++ ) {
		for ( int j = 0; j < S1.Len(p); j++ ) {
		  if ( pos >= exclude && pos <= N - exclude ) {
		    total_bases_this++;
		    total_cov_this += cov[pass][lib_id][p][j];
		  }
		  pos++;
		}
		if ( p < S1.Ntigs( ) - 1 ) pos += S1.Gap(p);
	      }
	      
              #pragma omp critical
	      {
		int64_t n_pairs = nfosmids_this[pass][lib_id];
		double cov = double(total_cov_this) / cov_frac;
		libdata[pass][lib_id].Add( n_pairs, total_bases_this, cov );
	      }
	    }
	  }
     }

     // Print table.

     size_t tot_reduced = 0;
     for( size_t ii=0; ii<scaffolds.size( ); ii++)
       tot_reduced += scaffolds[ii].ReducedLength( );

     out << Date( ) << ": counting reads in input" << endl;
     vec< vec<size_t> > n_input( 3 );   // frag, jump, and long_jump
     for (int pass=0; pass<3; pass++) {
       if ( ! lib_types[pass] ) continue;
       String s_type = pass == 0 ? "frag" : ( pass == 1 ? "jump" : "long_jump" );
       String orig_pairs_file = run_dir + "/" + s_type + "_reads_orig.pairs";
       PairsManager orig_pairs( orig_pairs_file );
       ForceAssertEq( orig_pairs.nLibraries( ), pairs[pass].nLibraries( ) );
       n_input[pass] = orig_pairs.getLibrarySizes( );   // number of pairs
       for (size_t ii=0; ii<n_input[pass].size( ); ii++)
	 n_input[pass][ii] *= 2;
     }     

     out << Date( ) << ": generating output table" << endl;
     vec< vec<String> > table;
     vec<String> line;
     line = MkVec( String( "type" ),
		   String( "lib_name" ),
		   String( "lib_stats" ),
		   String( "n_reads" ),
		   String( "%_used" ),
		   String( "scov" ),
		   String( "n_pairs" ),
		   String( "pcov" ) );
     table.push_back( line );

     vec<String> separator( 8, String( "" ) );
     table.push_back( separator );

     for (int pass=0; pass<3; pass++) {
       if ( ! lib_types[pass] ) continue;
       String s_type = pass == 0 ? "frag" : ( pass == 1 ? "jump" : "long_jump" );
       longlong nreads_tot = 0;
       longlong ninput_tot = 0;
       longlong tot_asslen = 0;
       longlong ntot_phys = 0;
       double ctot = 0.0;
       for (size_t lib_id=0; lib_id<pairs[pass].nLibraries( ); lib_id++) {
	 line.clear( );
	 String s_lname = pairs[pass].getLibraryName( lib_id );
	 int sep = pairs[pass].getLibrarySep( lib_id );
	 int dev = pairs[pass].getLibrarySD( lib_id );
	 String s_lstats = ToString( sep ) + " +/- " + ToString( dev );
	 String s_nreads = ToStringAddCommas( n_input[pass][lib_id] );
	 double pc_ass = double( assembled_count[pass][lib_id] ) / double( n_input[pass][lib_id] );
	 String s_pcreads = ToString( 100.0 * pc_ass, 1 );
	 double seq_cov = double( assembled_len[pass][lib_id] ) / double( tot_reduced );
	 String s_scov = ToString( seq_cov, 1 );
	 String s_npairs = ToStringAddCommas( libdata[pass][lib_id].NPairs( ) );
	 String s_pcov = ToString( libdata[pass][lib_id].PCov( ), 1 );
	 nreads_tot += assembled_count[pass][lib_id];
	 ninput_tot += n_input[pass][lib_id];
	 tot_asslen += assembled_len[pass][lib_id];
	 ntot_phys += libdata[pass][lib_id].NPairs( );
	 ctot += libdata[pass][lib_id].PCov( );

	 line = MkVec( s_type, s_lname, s_lstats, s_nreads, s_pcreads, s_scov, s_npairs, s_pcov );
	 table.push_back( line );
       }
       if ( pairs[pass].nLibraries( ) > 1 ) {   // no need for total if there is just one library
	 line.clear( );
	 String s_tot = "=== total ===";
	 String s_lib = "";
	 String s_nrtot = ToStringAddCommas( ninput_tot );
	 double tot_pc_ass = double( nreads_tot ) / double( ninput_tot );
	 String s_pcreads = ToString( 100.0 * tot_pc_ass, 1 );
	 String s_scov = ToString( double( tot_asslen ) / double( tot_reduced ), 1 );
	 String s_ntot_phys = ToStringAddCommas( ntot_phys );
	 String s_ctot = ToString( ctot, 1 );
	 
	 line = MkVec( s_type, s_tot, s_lib, s_nrtot, s_pcreads, s_scov, s_ntot_phys, s_ctot );
	 table.push_back( line );
       }
       table.push_back( separator );
     }

     if ( ! ARCHIVE ) PerfStat::log( ) << PerfStatBlockStart( "LibCoverage table" );
     if ( S != "all" ) out << "\nSelected supers: " << S << "\n";
     out << "\nLEGEND\n"
	 << "   n_reads:  number of reads in input\n"
	 << "   %_used:   % of reads assembled\n"
	 << "   scov:     sequence coverage\n"
	 << "   n_pairs:  number of valid pairs assembled\n"
	 << "   pcov:     physical coverage\n"
	 << endl;
     PrintTabular( out, table, 2, "llrrrrrr" );
     if ( ! ARCHIVE ) PerfStat::log( ) << PerfStatBlockStop( );

     out << Date( ) << ": done" << endl;
     if ( ARCHIVE ) cout << Date( ) << ": done" << endl;
     
}
