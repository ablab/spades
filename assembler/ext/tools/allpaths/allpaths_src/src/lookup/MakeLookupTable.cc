///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/// WARNING: This will not work correctly if there are blanks in the source
/// fasta files.  It would be better if we screened for this and aborted.

/// Build an index of all k-mers in a genome, as given by a collection of fasta
/// files.
/// \file MakeLookupTable.cc
///
/// MakeLookupTable.  Build an index of all k-mers in a genome, as given by a
/// collection of fasta files.
///
/// Note.  In the current design, the total number of fasta records in the files
/// should be small.
///
/// Arguments:
///
///      K: the k-mer size.  The maximum allowed value is 15.  Default: 12.
///
///      CHUNK_SIZE: the number of bases of genome which are indexed at one time.
///      Default: 100M.
///
///      CHUNK_OVERLAP: the overlap between chunks (because otherwise an alignment
///      across a chunk boundary might be lost).  Default: 1M.
///
///      SOURCE: an expression specifying the source fasta files for the genome
///      It may be enclosed in double quotes.  Files ending in ".gz" are treated
///      as gzipped.  (The expression is expanded by using the shell's "echo"
///      command.)  The number of bases in these files cannot exceed four billion.
///      Note: files ending in ".fastb" are assumed to be in fastb format.
///
///      OUT_HEAD: the name of the lookup table file prefix.  This code creates:
///      OUT_HEAD.lookup, OUT_HEAD.fastb, OUT_HEAD.fastamb.  The .lookup file is
///      the lookup table, the .fastb file contains the genome bases, and the
///      .fastamb file allows one to determine which bases were ambiguous.  For the
///      human genome, with default parameter choices, the .lookup file has size
///      about 14.5 GB.
///
///      IGNORE_SHORT: if set to True, suppress warning about source
///      sequences of length < K.  Otherwise, MakeLookupTable will
///      write a warning upon encountering such a sequence, which
///      cannot be hit by lookup table alignments, but is otherwise
///      harmless.
///
///      LOOKUP_ONLY: optional, default = False. If True, only lookup reference table
///      will be created (no fastb, fastamb files)
///
///      SORT_SOURCE: optional, default=True. If True, the source fasta file names
///      will be sorted before reading them in; otherwise they are read in order they
///      were specified on the command line.
///
///      CONTIGS: optional. Vector of (integer, zero-based) indexes of contigs to be
///      placed into the lookup table (and corresponding fastb, fastamb as well). The
///      indexes are with respect to the order, in which contigs will be actually
///      read in from the source fasta/fastb files (thus the result is affected by
///      the SORT_SOURCE option). If empty vector is specified (which is the
///      default), then all contigs are placed into the output files.

#include "Basevector.h"
#include "Bitvector.h"
#include "FastIfstream.h"
#include "FetchReads.h"
#include "FetchReadsAmb.h"
#include "FastaFileset.h"
#include "MainTools.h"
#include "lookup/LookupTable.h"
#include "lookup/LookupTableBuilder.h"
#include "feudal/IncrementalWriter.h"

int main( int argc, char *argv[] )
{
  RunTimeNoTraceback( );

  BeginCommandArguments;
  CommandArgument_UnsignedInt_OrDefault(K, 12);
  CommandArgument_UnsignedInt_OrDefault(CHUNK_SIZE, 100 * 1000 * 1000);
  CommandArgument_UnsignedInt_OrDefault(CHUNK_OVERLAP, 1000 * 1000);
  CommandArgument_String(SOURCE);
  CommandArgument_String_OrDefault_Doc(OUT_HEAD, SOURCE,
	   "Full qualified prefix for the output reference file(s)");
  CommandArgument_Bool_Abbr_OrDefault_Doc(LOOKUP_ONLY, LO, False,
		"If True, only reference lookup file will be created, without accompanying .fastb and .fastamb");
  CommandArgument_Bool_OrDefault_Doc(SORT_SOURCE, True,
	 "If True, multiple source files will be sorted before reading; otherwise they will be read in specified order");
  CommandArgument_Bool_OrDefault(IGNORE_SHORT, False);
  CommandArgument_IntSet_OrDefault_Doc(CONTIGS,"{}",
	    "Indices of contigs from source files to be added to the reference. Contigs are indexed globally, across all source files, in the order they will be read (see SORT_SOURCE)");
  CommandArgument_Bool_OrDefault(QUIET, False);
  CommandArgument_Bool_OrDefault_Doc(SIMPLE_SOURCE, False,
       "If True, SOURCE argument is not parsed into a file list.");
  EndCommandArguments;

  // Impose requirement on K.

  if ( K > 15 )
    FatalErr("The maximum allowed K value is 15.\n");

  // sort vector of contigs
  sort(CONTIGS.begin(),CONTIGS.end());

  // check that no stupid typos were made
  if ( CONTIGS.size() > 0 && Min(CONTIGS) < 0 ) {
    cerr << "Negative contig indexes are illegal. " << endl;
    exit(1);
  }


  // Get the genome fasta file names.

  vec<String> allfiles;
  if (SIMPLE_SOURCE) allfiles.push_back(SOURCE);
  else allfiles = AllFilesInSource(SOURCE);
  if (SORT_SOURCE) sort( allfiles.begin( ), allfiles.end( ), cmp_numeric );
  for ( int i = 0; i < (int) allfiles.size( ); i++ ) {
    if ( !QUIET ) cout << "using " << allfiles[i] << endl;
    if ( !IsRegularFile( allfiles[i] ) )
      FatalErr("Can't read file " << allfiles[i]);
    if ( !LOOKUP_ONLY && AreSameFile( OUT_HEAD + ".fastb", allfiles[i] ) )
      FatalErr("File OUT_HEAD.fastb, to be removed, is actually SOURCE file "
	       << allfiles[i]
	       << ".\nRemove the file by hand if that is your intent.");
  }

  // Generate .fastb and .fastamb files.

  int contig = 0; // global counter of contigs across all source files

  if ( !LOOKUP_ONLY ) {
    Remove( OUT_HEAD + ".fastb" ), Remove( OUT_HEAD + ".fastamb" );
    temp_file tmpfile( "/tmp/MakeLookupTable_tmp.XXXXXX" );
    vecbasevector v;
    vecbitvector b;
    String bFilename( OUT_HEAD + ".fastb" );
    String mbFilename( OUT_HEAD + ".fastamb" );
    IncrementalWriter<bitvector> mbWriter(mbFilename.c_str());
    IncrementalWriter<bvec> bWriter(bFilename.c_str());
    for ( int i = 0; i < (int) allfiles.size( ); i++ ) {
      v.clear( );
      b.clear( );
      if ( allfiles[i].Contains( ".gz", -1 ) ) {
	System( "gzip -c -d " + allfiles[i] + " > " + tmpfile );
	FastFetchReads( v, 0, tmpfile );
	FetchReadsAmb( b, tmpfile );
      } else if ( allfiles[i].Contains( ".fastb", -1 ) || IsGoodFeudalFile(allfiles[i]) ) {
	v.ReadAll( allfiles[i] );
	b.reserve( v.size() );
	for ( size_t j = 0; j < v.size( ); j++ ) {
	  static bitvector bi;
	  bi.clear().resize( v[j].size( ) );
	  b.push_back(bi);
	}
      } else {
	FastFetchReads( v, 0, allfiles[i] );
	FetchReadsAmb( b, allfiles[i] );
      }

      if ( CONTIGS.size() > 0 ) { // if we are selecting only a subset of contigs:
        // TODO: potentially dangerous truncation of index by to_remove
	vec<int> to_remove;
	for ( size_t j = 0 ; j < v.size() ; j++ ) { // for all contigs from current file
	  if ( -1 == BinPosition(CONTIGS,contig) ) {
	    // contig j from the current file was not requested; mark for discarding
	    to_remove.push_back(j);
	  }
	  contig++; // advance global contig counter
	}

	// remove contigs that were not asked for:
	v.RemoveByIndex(to_remove);
	b.RemoveByIndex(to_remove);
      }

      bWriter.add(v.begin(),v.end());
      mbWriter.add(b.begin(),b.end());
    }

    // now we know what the total number of contigs in source files was.
    // sanity check:
    if ( CONTIGS.size() > 0 && Max(CONTIGS) >= contig ) {
      cerr << "Illegal contig index was specified. Requested: " << Max(CONTIGS) <<
	"; total number of contigs available in source files: " <<
	contig << endl;
      exit(1);
    }

    bWriter.close();
    mbWriter.close();
  }

  // Set up to create .lookup file.

  Remove( OUT_HEAD + ".lookup" );
  int fd = Open( OUT_HEAD + ".lookup", O_WRONLY | O_CREAT );

  // Count number of bases, possibly slightly overcounting.
  // NOTE: this is the total count for all contigs found in
  // source files; currently does not take into account that only a subset of contigs
  // could have been requested! It's a HUGE overcount.
  const ulonglong max_bases(4000000000u);
  ulonglong base_count = 0;
  ulonglong contig_count = 0;
  for ( int i = 0; i < allfiles.isize( ); i++ ) {
    if ( allfiles[i].Contains( ".fastb", -1 ) || IsGoodFeudalFile(allfiles[i]) ) {
      vecbasevector ref(allfiles[i]);
      contig_count += ref.size();
      for ( unsigned int j = 0 ; j < static_cast<unsigned int>(ref.size()) ; j++ ) {
	base_count+=ref[j].size();
      }

      // The old code below fails to count bases correctly, because of the memory alignment
      // issues: basevector uses 2 bits per base but allocates memory by unsigned int (4-byte)
      // units. MastervecFileRawCount returns the minimum number of such allocation units
      // needed to store the actual bases, and it is the total number of bases that can be stored
      // in these units that is computed belo. Hence the base count for each contig can be off
      // by as much as 15 bases  ( (4 bytes *8 bits )/2 bits_per_base - 1 ).
      /*
      // basevectors are Feudal vectors of unsigned int, packed with 4 bases per byte
         base_count += MastervecFileRawCount(allfiles[i]) * sizeof(unsigned int) * 4;
         contig_count += MastervecFileObjectCount(allfiles[i]);
      */
    } else {
      // Cat or zcat file so we can deal equivalently with gzipped files
      String command = "cat ";
      if ( allfiles[i].Contains( ".gz", -1 ) ) command = "z" + command;
      command += allfiles[i];
      // And count the printing characters on the non-label lines
      command += " | grep -v '^>' | tr -d -c '[:print:]' | wc -c";
      String count = StringOfOutput(command);
      base_count += count.Int();
    }
  }
  if ( contig_count > 0 ) cout << "Total " << contig_count << " contigs in the reference file(s)" << endl;
  cout << "Total " << base_count << " bases in the reference sequence(s)" << endl;
  flush(cout);

  if ( base_count >= max_bases )
    FatalErr("Your input files contain more than four billion bases.  Abort.");

  // Set up lookup table.

  lookup_table look(fd);
  look.SetK(K);
  look.SetChunkParams( CHUNK_SIZE, CHUNK_OVERLAP );

  // Read in the bases.
  vec<char> bases(base_count);
  base_count = 0;
  int fastb = 0;
  for ( int i = 0; i < (int) allfiles.size( ); i++ )
    if ( allfiles[i].Contains( ".fastb", -1 ) ) ++fastb;
  if ( fastb > 0 && fastb < (int) allfiles.size( ) )
    FatalErr("All or none of the source files must be in fastb format. Abort.");

  contig = 0; // reset contig counter

  for ( int i = 0; i < (int) allfiles.size( ); i++ ) {
    if ( allfiles[i].Contains( ".fastb", -1 ) || IsGoodFeudalFile(allfiles[i]) ) {
      vecbasevector m;
      m.ReadAll( allfiles[i] );

      // TODO: potentially dangerous truncation of index by records
      vec<int> records; // will keep original index of the contigs (for annotation)
      if ( CONTIGS.size() > 0 ) { // if only a subset of contigs was requested:
	vec<int> to_remove;
	for ( size_t j = 0 ; j < m.size(); j++ ) { // for each contig from current file:
	  if ( -1 == BinPosition(CONTIGS,contig) ) {
	    // check absolute index of the contig (across all source files).
	    // if contig j from the file was not requested, mark it for discarding:
	    to_remove.push_back(j);
	  } else {
	    records.push_back(j); // store the contig's index
	  }
	  // advance global contig counter; note: if CONTIGS.size()==0, we don't
	  // need this counter at all, so it's ok to increment it inside if() {}:
	  contig++;
	}
	m.RemoveByIndex(to_remove); // keep only requested contigs
      } else {
	for ( size_t j = 0 ; j < m.size() ; j++ ) records.push_back(j); // all records
      }

      //      int record = 0;
      for ( size_t j = 0; j < m.size( ); j++ ) {
	unsigned int n = m[j].size( );
	// here we use 'records' initialized above: alternative contig name
	// is <source_filename>:<index_in_source_file>:
	look.AddContigName(
			   "contig_" + ToString(j), allfiles[i], records[j] );
	look.AddContigStart(base_count);
	//	++record;
	look.AddContigSize(n);
	if ( n < K  && !IGNORE_SHORT ) {
	  if ( !QUIET ) cout << "Warning: Contig \"" << look.LastContigName( )
	       << "\"\ncannot be aligned to because its length is " << n
	       << ", which is less than K.\n";
	} else if ( n > 2000000000u ) {
	  FatalErr("Size of contig \""
		   << look.LastContigName( ) << "\" is " << n
		   << ", which is too large "
		   << "(max value = 2,000,000,000).\n"); // FATAL ERROR!
	}
	// add bases from the current contig to the concatenated seq. vector
	for ( unsigned int k = 0; k < n; k++ )
	  bases[ base_count++ ] = as_base( m[j][k] );
      }
      continue; // done with current fastb file, go get next
    }
    // if we got here, the current file is (should be) .fasta
    String line, command = "cat " + allfiles[i];
    if ( allfiles[i].Contains( ".gz", -1 ) ) command = "z" + command;
    fast_pipe_ifstream in(command);
    //    Bool first = True;

    // are we currently reading a contig? this flag is always true when all
    // contigs are being read; it will change states when a subset of contigs
    // is read and some contigs need to be skipped (then reading=false happens):
    Bool reading = False;
    int record = 0; // counts original records in the source file

    // preprocessing (moved out from the loop below; this code was executed
    // only once anyway, but required additional flags to ensure single execution
    // while inside the loop):

    getline( in, line );  // read the first line
    if ( in.fail() ) {
      FatalErr(allfiles[i] << "is empty. Abort.");
    }
    if ( !line.Contains( ">", 0 ) ) { // naive check. is this fasta?
      FatalErr(allfiles[i] << " is not in fasta format.  Abort.");
    }

    if ( CONTIGS.size() == 0 || -1 != BinPosition(CONTIGS,contig) ) {
      // if all contigs were requested, or if the current contig (first
      // contig in the current source file) was explicitly requested:
      look.AddContigName( line.After( ">" ), allfiles[i], record );
      look.AddContigStart(base_count);
      reading = True;
    } else {
      // we need to skip this contig:
      reading = False;
    }

    while(1) {
      getline( in, line );  // read the next line

      if ( line.Contains( ">", 0 ) )  { // new record (contig) is found

	if ( reading )  {
	  // if we were reading, not skipping, the contig that just ended:
	  unsigned int cs = base_count - look.LastContigStart( );

	  if ( cs < K && !IGNORE_SHORT ) {
	    if ( !QUIET ) cout << "Warning: Contig \"" << look.LastContigName( )
		 << "\"\ncannot be aligned to because its length is " << cs
		 << ", which is less than K.\n";
	  } else if ( cs > 2000000000u ) {
	    FatalErr("Size of contig \""
		     << look.LastContigName( ) << "\" is " << cs
		     << ", which is too large "
		     << "(max value = 2,000,000,000).\n"); // FATAL ERROR!
	  }
	  look.AddContigSize(cs); // the name of the contig that just ended was set
	                          // when it started (from its ">" line), now set size
	}
        // a new contig is about to start; increment global contig counter and
	// local record counter:
	contig++;
	record++;
	if ( 0 == CONTIGS.size() || -1 != BinPosition(CONTIGS,contig) ) {
	  // if we need to read all contigs, or if the contig that has just
	  // started was explicitly requested, add name/start and mark for reading:
	  look.AddContigName( line.After( ">" ), allfiles[i], record );
	  look.AddContigStart(base_count);
	  reading = True;
	} else {
	  reading = False; // the contig we joust found has to be skipped
	}
      }
      else {
	// we are inside the record; if we are not skipping the contig, read bases:
	if ( reading ) {
	  for ( int j = 0; j < (int) line.size( ); j++ ) {
	    if (isprint(line[j])) // ignore nonprinting characters
	      bases[ base_count++ ] = line[j];
	  }
	}
      }
      if ( in.fail( ) ) break;
    }

    // the file has ended; let's add last contig size, if that contig was requested:
    if ( reading ) {
      look.AddContigSize( base_count - look.LastContigStart( ) );
    }
  } // end of 'for ( allfiles )' loop across source file names


  // The original base count might have been a slight overestimate, so fix it
  ForceAssertLe(base_count, bases.size());
  bases.resize(base_count);

  // Build the table
  BuildTableFromContigs(look, bases, K, CHUNK_SIZE, CHUNK_OVERLAP);

  cout << Date( ) << ": Done with MakeLookupTable!" << endl;
  return 0;
}
