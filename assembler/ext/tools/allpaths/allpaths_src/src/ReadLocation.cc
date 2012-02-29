///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include <map>
#include <numeric>

#include "CoreTools.h"
#include "ReadLocation.h"
#include "ReadPairing.h"
#include "VecString.h"
#include "VecTemplate.h"
#include "system/file/FileReader.h"

void read_location::ForceInBounds( ostream * out_ptr )
{    if ( (int) LengthOfRead( ) > LengthOfContig( ) )
     {    if ( out_ptr )
          {    *out_ptr << "Warning! ForceInBounds: read length = " 
                    << LengthOfRead( ) << ", longer than contig length = " 
                    << LengthOfContig( ) << ", for read = " << ReadId( ) 
                    << ", contig = " << Contig( ) << "\n";    }
          return;    }
     if ( StartOnContig( ) < 0 )
     {    if ( out_ptr && StartOnContig( ) < -10 )
               *out_ptr << "Warning: start location for read " << ReadId( ) 
                       << " on contig " << Contig( )
                       << " would have been " << StartOnContig( ) << ".\n";
          SetStartOnContig(0);    }
     if ( StopOnContig( ) > LengthOfContig( ) )
     {    int extra = StopOnContig( ) - LengthOfContig( );
          if ( out_ptr && extra > 10 )
               *out_ptr << "Warning: stop location for read " << ReadId( )
                        << " on contig " << Contig( )
                        << " would have been " << extra 
                        << " beyond the the contig end.\n";
          SetStartOnContig( StartOnContig( ) - extra );
     }    
}

istream& operator>>( istream& s, vec<read_location>& v )
{    int n;
     s >> n;
     v.resize(n);
     char c;
     s.get(c);
     for ( int i = 0; i < n; i++ )
          s.read( (char*) &v[i], sizeof(read_location) );
     return s;    }

void ReadsToLocations( const String& locs_file, const vec<int>& ids,
     vec< vec<read_location> >& locs )
{    String nlocs;
     {    Ifstream( in, locs_file );
          in >> nlocs;    }
     int head = nlocs.size( ) + 1;
     FileReader fr( locs_file.c_str() );
     FileReader fri( (locs_file + "_indexr").c_str() );
     locs.resize( ids.size( ) );
     int four = sizeof(int);
     for ( int i = 0; i < ids.isize( ); i++ )
     {    int id = ids[i];
          if ( id < 0 ) 
          {    locs[i].clear( );
               continue;    }
          int idx;
          fri.seek( four + id * four );
          fri.read( &idx, four );
          if ( idx >= 0 ) 
          {    locs[i].resize(1);
               fr.seek( head + idx * sizeof(read_location) );
               fr.read( &locs[i][0], sizeof(read_location) );    }
          else if ( idx < -1 )
          {    fri.seek( four + (-idx * four) );
               int nplacements;
               fri.read( &nplacements, sizeof(int) );
               vec<int> placements(nplacements);
               fri.read( &placements[0], nplacements * sizeof(int) );
               locs[i].resize(nplacements);
               for ( int j = 0; j < nplacements; j++ )
               {    fr.seek( head + placements[j] * sizeof(read_location) );
                    fr.read( &locs[i][j], sizeof(read_location) );    }    }
          else locs[i].clear( );    }    }

void WriteLocs( const String& locs_file , const vec<read_location>& v, 
     const int num_contigs, int nreads, Bool write_locs )
{    
     // Check for bad data.

     for ( int i = 0; i < v.isize( ); i++ )
          ForceAssertGe( v[i].Contig( ), 0 );

     // Write read locations to locs_file.

     if (write_locs)
     {    ofstream wr_strm(locs_file.c_str(), ios::out|ios::binary);
          wr_strm << v.size( ) << "\n";
          for ( unsigned int i = 0; i < v.size( ); i++ )
               wr_strm.write( (char*) &v[i], sizeof(read_location) );    }

     // Generate index by reads.  This is a vector u of ints, as follows: first 
     // there are nreads ints.  For a given read, if the read is unplaced, the entry 
     // is -1.  If the read is placed uniquely, the entry is the index of its read 
     // location.  If the read is multiply placed, then the entry is -x, where x 
     // is the index of a multiple read placement block in u (to be described).  
     // After the nreads ints in u come additional entries, for the multiply placed 
     // reads.  These additional entries are in blocks.  The block starts with the 
     // number of multiple placements (of a given read), then has the indices of its 
     // read locations.

     String locs_index_filer = locs_file + "_indexr";
     if ( nreads < 0 ) Remove(locs_index_filer);
     else
     {
          // First count how many times each read is placed.  We need this to make
          // the index the right size.  Then build the index.

          vec<int> nplaced( nreads, 0 );
          int nlocs = v.size( );
          for ( int i = 0; i < nlocs; i++ )
               ++nplaced[ v[i].ReadId( ) ];
          int total_entries = nreads;
          for ( int i = 0; i < nreads; i++ )
               if ( nplaced[i] > 1 ) total_entries += (1 + nplaced[i]);
          vec<int> indexr( total_entries, -1 );
          int ptr = nreads;
          for ( int i = 0; i < nreads; i++ )
          {    if ( nplaced[i] > 1 )
               {    indexr[i] = -ptr;
                    indexr[ptr] = nplaced[i];
                    ptr += (1 + nplaced[i]);    }    }
          for ( int i = 0; i < nlocs; i++ )
          {    int id = v[i].ReadId( );
               if ( nplaced[id] == 1 ) indexr[id] = i;
               else if ( nplaced[id] > 1 )
               {    for ( int j = -indexr[id] + 1; ; j++ )
                    {    if ( indexr[j] < 0 ) 
                         {    indexr[j] = i;
                              break;    }    }    }    }
          BinaryWrite( locs_index_filer, indexr );    }

     // Generate index by contigs.

     String locs_index_file = locs_file + "_index";
     if ( num_contigs < 0 ) Remove(locs_index_file);
     else
     {
          if ( v.empty() ) 
               FatalErr( "Empty read_location vector passed to WriteLocs." );

          // determine the index boundaries between contigs 
          // loc_index stores the beginning read index for each contig
          // NOTE: the loc_index array index is the contig number

          // Initialize all boundaries to -1.
          vec<int>  by_contig_index(num_contigs+1,-1);

          // Walk backwards through the locs, updating the index for
          // each contig.  The entry by_contig_index[i] will end up
          // containing the index of the first loc for contig i
          // (presuming contig i has reads in it).

          int lastContig = v.back().Contig();
          if ( num_contigs <= lastContig )
          {
            FatalErr( "Incorrect num_contigs passed to WriteLocs: " <<
                      "num_contigs=" << num_contigs <<
                      " and v.back().Contig()=" << lastContig << "." );
          }

          for ( int i = v.size()-1; i >= 0; --i )
          {
            int thisContig = v[i].Contig();
            if ( thisContig > lastContig )
            {
              FatalErr( "Unsorted read_location vector passed to WriteLocs "
                        << "with request for index by contigs." );
            }
            by_contig_index[ thisContig ] = i;
            lastContig = thisContig;
          }

          // Set the start of empty contigs to be that of the start of
          // the next non-empty contig.  Start with last as the size
          // of the locs vector so the last index is correct.
          int last = v.size();
          for ( int i = by_contig_index.size()-1; i >= 0; --i )
            if ( by_contig_index[i] < 0 )
              by_contig_index[i] = last;
            else
              last = by_contig_index[i];
       
          // write by_contig_index to the file locs_file_index
     
          ofstream ind_wr_strm(locs_index_file.c_str(), ios::out|ios::binary);
          for ( unsigned int i = 0; i < by_contig_index.size(); i++ )
	      ind_wr_strm.write((char *) &by_contig_index[i], sizeof(unsigned int));
          ind_wr_strm.close();    }
}

// Helper routine to read indices.
void ReadLocIndex( const String & locs_file , vec<int> & start_loc_of_contig )
{
  start_loc_of_contig.clear();
  
  if ( !IsRegularFile(locs_file) )
    FatalErr( locs_file << " is not a regular file." );

  // the index file
  String locs_index_file(locs_file);
  locs_index_file += "_index";
  if ( !IsRegularFile(locs_index_file) )
    FatalErr( locs_index_file << " is not a regular file." );

  // read locs index file and put in loc_index
  //
  // open file and go to end, easy way of getting size
  ifstream index_rd_strm(locs_index_file.c_str(), ios::in|ios::binary|ios::ate );
  if ( !index_rd_strm )                          
    FatalErr( "Problem reading " << locs_index_file << "." );    
  streampos size = index_rd_strm.tellg(); 
  ForceAssert(size > 0 );
  long num_loc_index = size/sizeof(int);    
  index_rd_strm.seekg(0, ios::beg);

  start_loc_of_contig.resize( num_loc_index );
  index_rd_strm.read((char *) &start_loc_of_contig[0], size);
  index_rd_strm.close();
}


void ReadLocs( const String & locs_file , const vec<int> & contig_ids, vec<read_location>& v ) 
{
     v.clear();

     if ( !IsRegularFile(locs_file) )
       FatalErr( locs_file << " is not a regular file." );

     if ( contig_ids.size() == 0 ) return;

     vec<int> loc_index;
     ReadLocIndex( locs_file, loc_index );

     // open locs file, jump past first integer (size)
     ifstream rd_strm(locs_file.c_str(), ios::in|ios::binary );
     if ( !rd_strm )                          
         FatalErr( "Problem reading " << locs_file << "." ); 
     unsigned int n;
     rd_strm >> n; // bypass size
     char c;
     rd_strm.get(c); // bypass EOL
     streampos start_of_data = rd_strm.tellg();

     // put only locs with contig_ids in v
     unsigned int tot_reads = 0;
     unsigned int read_loc_size = sizeof(read_location);
     for ( unsigned int i = 0; i < contig_ids.size(); i++ ){ 
       int contig = contig_ids[i];
       ForceAssertLt( contig + 1, (int) loc_index.size() );
       ForceAssertLt( loc_index[contig], loc_index[contig+1] );
       streampos pos = loc_index[contig]*read_loc_size + start_of_data;
       unsigned int num_reads = loc_index[contig+1] - loc_index[contig];
       unsigned int len = num_reads*read_loc_size;
       tot_reads += num_reads;
       rd_strm.seekg(pos, ios::beg);
       unsigned int prev_size = v.size();
       v.resize(tot_reads);
       rd_strm.read((char *) &v[prev_size], len);
       ForceAssertEq( v[prev_size].Contig(), contig ); 
     }
     rd_strm.close();
}

void NumLocs( const String& locs_file, vec<int> &num_locs_per_contig )
{
     if ( !IsRegularFile(locs_file) )
       FatalErr( locs_file << " is not a regular file." );
     
     ReadLocIndex( locs_file, num_locs_per_contig );
     
     adjacent_difference( num_locs_per_contig.begin(), 
                          num_locs_per_contig.end(),
                          num_locs_per_contig.begin() );

     if ( ! num_locs_per_contig.empty() )
     {
       copy( num_locs_per_contig.begin() + 1,
             num_locs_per_contig.end(),
             num_locs_per_contig.begin() );
       num_locs_per_contig.resize( num_locs_per_contig.size() - 1 );
     }

     for ( unsigned int contig_id = 0; contig_id < num_locs_per_contig.size(); ++contig_id )
       ForceAssertGe( num_locs_per_contig[ contig_id ], 0 );
}
  

void AnnotateReadLocations( vec<read_location>& data, String run_dir,
     Bool orig, String human_file, const vec<Bool>& contigs_to_use, Bool gzip,
     Bool untrimmed, Bool show_divider )
{    
     Remove( human_file ), Remove( human_file + ".gz" );

     String run_dirx = ( !orig ? run_dir : run_dir + "/orig" );

     int N = MastervecFileObjectCount( run_dir + "/reads.fastb" );

     vec<Bool> is_transposon;
     if (orig) READX( run_dir + "/reads.is_transposon", is_transposon );

     vec< vec<int> > datax;
     datax.resize(N);
     for ( unsigned int i = 0; i < data.size( ); i++ )
     {    ForceAssertLt( data[i].ReadId( ), N );
          datax[ data[i].ReadId( ) ].push_back(i);    }

     vecString ids( run_dirx + "/reads.ids" );

     vec<int> left_trim, right_trim;
     if (untrimmed)
     {    Ifstream( lr_trim, run_dir + "/reads.trim_lr" );
          left_trim.reserve( N );
          right_trim.reserve( N );
          while(1)
          {    int n;
               lr_trim >> n;
               if ( !lr_trim ) break;
               left_trim.push_back(n);
               lr_trim >> n;
               right_trim.push_back(n);    }    }

     vec<read_pairing> pairs;
     ReadPairsFile( run_dirx + "/reads.pairto", pairs );
     vec<int> pairs_index(N, -1);
     for ( unsigned int i = 0; i < pairs.size( ); i++ )
          if ( pairs[i].Alive( ) )
               pairs_index[ pairs[i].id1 ] = pairs_index[ pairs[i].id2 ] = i;

     ostream* s;
     procbuf* pb;
     if ( gzip ) {
       String command = "gzip -1 > " + human_file + ".gz"; 
       pb = new procbuf( command.c_str(), ios::out );
       s = new ostream( pb );
     }
     else {
       pb = 0;
       s = new ofstream(human_file.c_str(), ios::out|ios::binary);
     }

     ostream& human_locs = *s;
     
     vec<read_location>::iterator rl_set_begin = data.begin();
     vec<read_location>::iterator rl_set_end = data.begin();

     Bool first_contig = True;
     while ( ++rl_set_end <= data.end() &&
	     rl_set_begin < data.end() )
     {    
       if ( rl_set_end != data.end() &&
	    rl_set_end->Contig() == rl_set_begin->Contig() )
	 continue;

       if ( contigs_to_use.size( ) > 0 && !contigs_to_use[ rl_set_begin->Contig( ) ] )
       {
	 rl_set_begin = rl_set_end;
	 continue;
       }

       vec<read_location> rl_set( rl_set_end - rl_set_begin );
       copy( rl_set_begin, rl_set_end, rl_set.begin() );

       if ( untrimmed )
       {
	 for ( unsigned int i = 0; i < rl_set.size(); ++i )
	 {
	   read_location& rl = rl_set[i];
	   int id = rl.ReadId( );

	   rl.SetLengthOfRead( rl.LengthOfRead() + left_trim[id] + right_trim[id] );
	   if ( rl.OrientationOnContig() == ForwardOr )
	     rl.SetStartOnContig( rl.StartOnContig() - left_trim[id] );
	   else
	     rl.SetStartOnContig( rl.StartOnContig() - right_trim[id] );
	 }
       }

       Sort(rl_set);

       if ( show_divider && !first_contig )
       {    human_locs << "------------------------------------------"
                 << "------------------------------------------\n\n";    }
       first_contig = False;

       for ( unsigned int i = 0; i < rl_set.size(); ++i )
       {
          read_location& rl = rl_set[i];
          int id = rl.ReadId( );
          human_locs << ids[id] << "\n" << rl;

          int id_partner = pairs_index[id] >= 0 ?
	    pairs[ pairs_index[id] ].Partner(id) : -1;

          if ( id_partner == -1 ) human_locs << "(not paired)\n\n";    
          else 
          {    
	    human_locs << "(paired to " << id_partner << " = " << ids[id_partner];
	    
	    const vec<int>& pl = datax[id_partner];
	    
	    if ( pl.size( ) == 0 ) human_locs << ", not in the assembly)\n\n";

	    else if ( pl.size( ) > 1 ) human_locs << ", multiply located)\n\n";
	    
	    else
	    {
	      read_location& partner_in_data = data[ pl[0] ];
	      if ( partner_in_data.Contig( ) != rl.Contig( ) )
		human_locs << ", in contig " << partner_in_data.Contig() << ")\n\n";
	      else
	      {
		unsigned int rlp_set_index = 0; 
		for ( ; rlp_set_index < rl_set.size(); ++rlp_set_index )
		  if ( rl_set[rlp_set_index].ReadId() == id_partner )
		    break;
		read_location& rlp = rl_set[rlp_set_index];

		if ( rlp_set_index < i )
		  human_locs << ", above in this contig)\n\n";

		else if ( ( !orig || !is_transposon[ rl.ReadId( ) ] ) &&
			  rl.OrientationOnContig( ) == ForwardOr &&
			  rlp.OrientationOnContig( ) == ReverseOr )
		{ 
		  int insert_size = rlp.StopOnContig( ) - rl.StartOnContig( ) + 1;
                  int expected = pairs[ pairs_index[id] ].sep;
                  int diff = rlp.StartOnContig( ) - rl.StopOnContig( ) - expected;
                  int sd = pairs[ pairs_index[id] ].sd;
                  ForceAssert( sd != 0 );
                  float discrep = float(diff) / float(sd);
		  human_locs << ", below in this contig,\n" << "insert size ~ " 
                       << insert_size << " = expected";
                  if ( diff >= 0 ) 
                       human_locs << " + " << diff << " [+" << setprecision(3)
                            << discrep << " sd]";
                  else human_locs << " - " << -diff << " [-" << setprecision(3)
                            << -discrep << " sd]";
                  human_locs << ")\n\n";
		}

		else if ( ( orig && is_transposon[ rl.ReadId() ] ) &&
			  rl.OrientationOnContig( ) == ReverseOr &&
			  rlp.OrientationOnContig( ) == ForwardOr )
		{ 
		  int sep = rlp.StartOnContig( ) - rl.StopOnContig( ) + 1;
		  if ( untrimmed )
		    sep += left_trim[rl.ReadId()] + left_trim[rlp.ReadId()];
		  human_locs << ", below in this contig, " 
			     << "separation after trimming ~ " 
			     << sep << ")\n\n";
		}
		else 
		  human_locs << ", below in this contig, illogically)\n\n"; 
	      }
	    }
	  }
       }


       rl_set_begin = rl_set_end;
     }
     
     delete s;
     delete pb;
}

void AnnotateReadLocations2( const String& run_dir, const String& sub_dir,
     const vec<int>& tig_ids, ostream& out )
{    String run_dirx = run_dir + "/orig";
     vec<read_location> locs;
     String locs_file = sub_dir + "/mergedcontigs_orig.locs";
     ReadLocs( locs_file, tig_ids, locs );
     int n = locs.size( );
     vec<int> read_ids(n), partner_ids(n, -1), all_ids;
     map<int, int> to_ids_plus;
     for ( int i = 0; i < n; i++ )
     {    read_ids[i] = locs[i].ReadId( );
          to_ids_plus[ read_ids[i] ] = i + 1;    }
     vec<read_pairing> pairs;
     ReadsToPairs( run_dirx, read_ids, pairs );
     for ( int i = 0; i < n; i++ )
          if ( pairs[i].Alive( ) ) partner_ids[i] = pairs[i].Partner( read_ids[i] );
     vec< vec<read_location> > partner_locs;
     ReadsToLocations( locs_file, partner_ids, partner_locs );
     all_ids.reserve( 2*n ), all_ids.append(read_ids);
     for ( int i = 0; i < n; i++ )
          if ( partner_ids[i] >= 0 ) all_ids.push_back( partner_ids[i] );
     UniqueSort(all_ids);
     map<int, int> to_all;
     for ( int i = 0; i < all_ids.isize( ); i++ )
          to_all[ all_ids[i] ] = i;
     vec<Bool> is_transposon;
     BinaryReadSubset0( run_dir + "/reads.is_transposon", all_ids, is_transposon );
     vecString readnames;
     readnames.Read( run_dirx + "/reads.ids", all_ids );
     int last_contig = -1;
     for ( int i = 0; i < n; i++ )
     {    if ( last_contig >= 0 && locs[i].Contig( ) != last_contig )
          {    out << "------------------------------------------"
                    << "------------------------------------------\n\n";    }
          last_contig = locs[i].Contig( );
          const read_location& rl = locs[i];
          int id = rl.ReadId( );
          out << readnames[ to_all[id] ] << "\n" << rl;
          int id_partner = partner_ids[i];
          if ( id_partner == -1 ) out << "(not paired)\n\n";    
          else 
          {    out << "(paired to " << id_partner << " = " 
                    << readnames[ to_all[id_partner] ];
	       const vec<read_location>& pl = partner_locs[i];
	       if ( pl.empty( ) ) out << ", not in the assembly)\n\n";
	       if ( pl.size( ) > 1 ) out << ", multiply located)\n\n";
               if ( pl.size( ) != 1 ) continue;
	       const read_location& rlp = pl[0];
	       int k = to_ids_plus[id_partner] - 1;
	       if ( rlp.Contig( ) != rl.Contig( ) )
	            out << ", in contig " << rlp.Contig( ) << ")\n\n";
               else if ( k < i ) out << ", above in this contig)\n\n";
	       else if ( !is_transposon[i] && rl.Fw( ) && rlp.Rc( ) )
	       {    int insert_size = rlp.StopOnContig( ) - rl.StartOnContig( ) + 1;
                    int expected = pairs[i].sep;
                    int diff = rlp.StartOnContig( ) - rl.StopOnContig( ) - expected;
                    int sd = pairs[i].sd;
                    ForceAssert( sd != 0 );
                    float discrep = float(diff) / float(sd);
		    out << ", below in this contig,\n" << "insert size ~ " 
                         << insert_size << " = expected";
                    if ( diff >= 0 ) 
                         out << " + " << diff << " [+" << setprecision(3)
                              << discrep << " sd]";
                    else out << " - " << -diff << " [-" << setprecision(3)
                              << -discrep << " sd]";
                    out << ")\n\n";    }
   	       else if ( ( is_transposon[i] ) && rl.Rc( ) && rlp.Fw( ) )
	       {    int sep = rlp.StartOnContig( ) - rl.StopOnContig( ) + 1;
		    out << ", below in this contig, " 
		         << "separation after trimming ~ " << sep << ")\n\n";    }
    	       else 
               {    out << ", below in this contig, illogically)"
                         << "\n\n";     }    }    }    }

bool operator==(const read_location& oldRead, const read_location& newRead) {
       return ((oldRead.ReadId() == newRead.ReadId())
	       && (oldRead.LengthOfRead() == newRead.LengthOfRead())
	       && (oldRead.Contig() == newRead.Contig())
	       && (oldRead.StartOnContig() == newRead.StartOnContig())
	       && (oldRead.StopOnContig() == newRead.StopOnContig())
	       && (oldRead.OrientationOnContig() == newRead.OrientationOnContig())
	       && (oldRead.LengthOfContig() == newRead.LengthOfContig())
	       && (oldRead.OnSupers() == newRead.OnSupers())
	       && (oldRead.SuperContig() == newRead.SuperContig())
	       && (oldRead.StartOnSuperContig() == newRead.StartOnSuperContig())
	       && (oldRead.StopOnSuperContig() == newRead.StopOnSuperContig())
	       && (oldRead.LengthOfSuperContig()
		   == newRead.LengthOfSuperContig()));   }

ostream& operator<<(ostream &o, const read_location &r)
{    o << "rd " << r.ReadId( ) << " (len=" << r.LengthOfRead( ) << 
          ") --> " << (r.OrientationOnContig( ) == ForwardOr ? "" : "rc ")
          << "tig " << r.Contig( ) << " (len=" << r.LengthOfContig( )
          << ") @ " << r.StartOnContig( );
     if ( r.OnSupers( ) )
          o << " --> ^" << r.SuperContig( ) << " (len=" << r.LengthOfSuperContig( ) 
               << ") @ " << r.StartOnSuperContig( );
     return o << "\n";    }

bool operator%(const read_location& oldRead, const read_location& newRead) {
       if ((oldRead.ReadId() == newRead.ReadId())
	       && (oldRead.LengthOfRead() == newRead.LengthOfRead())
	       && (oldRead.Contig() == newRead.Contig())
	       && (oldRead.StartOnContig() == newRead.StartOnContig())
	       && (oldRead.StopOnContig() == newRead.StopOnContig())
	       && (oldRead.OrientationOnContig() == newRead.OrientationOnContig())
	       && (oldRead.LengthOfContig() == newRead.LengthOfContig()))
	   return 0;
       else return 1;   }

BINARY2_DEF(read_location_short);

