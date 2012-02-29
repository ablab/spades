///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*******************************************************************************
 *
 * Class: PairsManager
 *
 *
 * A PairsManager is a collection of all of the read-read pairings in a dataset.
 *
 * Each read pair is described fully by:
 * -- Two read IDs.  These are assumed to correspond to the indices in a
 *    vecbasevector or vecKmerPath or whatever.
 * -- A library ID.  This number is indexed into an internal set of libraries.
 *
 * Unpaired reads are not listed; hence the number of pairs may be smaller than
 * half the number of reads.
 *
 * Each library is contained in a PM_Library struct, which contains information
 * about the libraries: separation and standard deviation, plus a library name.
 *
 * The PairsManager class is designed to replace the older read_pairing class
 * (ReadPairing.h) throughout the RunAllPathsLG pipeline.  A PairsManager uses
 * less memory and smaller file sizes than a vec<read_pairing>, and can handle
 * larger datasets (more than 2^31 reads).
 *
 * We have almost completely succeeded in replaing all read_pairings with
 * PairsManagers throughout the RunAllPathsLG pipeline.  This class contains
 * some conversion functions to simplify the process of drop-in replacement.
 *
 *
 *
 * MEMORY USAGE AND FILE SIZE:
 * 17 bytes per pair, plus a small constant amount.
 *
 * NUMERIC LIMITATIONS:
 * Number of reads/pairs must be less than max(longlong) = 2^63 =~ 1e19
 * Number of libraries must be less than max(unsigned char) = 127
 * Read separation and stdev must be less than max(int) = 2^31; can be negative
 *
 *
 *
 * Josh Burton
 * August 2009
 *
 ******************************************************************************/



#ifndef PAIRS_MANAGER_H
#define PAIRS_MANAGER_H

#include "ReadPairing.h" // read_pairing
#include "String.h"
#include "Vec.h"
#include "feudal/BinaryStream.h"
#include <cstddef>


class PairsManager; // forward declaration


/******************************************************************************
 *
 *            PM_Library: Library info, for use with PairsManager
 *
 ******************************************************************************/
struct PM_Library {
  
  const PairsManager * _pm; // Pointer to the PairsManager containing this lib
  String _name; // Should be a unique identifier, e.g., "Solexa-14967"
  int _sep, _sd; // Pair separation and standard deviation
  
  
  // PM_Library constructor.  Provides a default name if none is given.
  PM_Library( const PairsManager * pm, const int sep, const int sd, const String & name = "" );
  
  
  friend bool operator==( const PM_Library & lib1, const PM_Library& lib2 );
  friend bool operator!=( const PM_Library & lib1, const PM_Library& lib2 );
};


/******************************************************************************
 *
 *          PM_LibraryStats: used to export library information from PM class  
 *
 ******************************************************************************/

struct PM_LibraryStats {
  PM_LibraryStats() {
    name = "";
    sep = sd = 0;
    min_len = max_len = mean_len = 0;
    n_bases = n_reads = 0;
  }

  String name;
  int32_t sep, sd;
  uint32_t min_len, max_len, mean_len;
  uint64_t n_bases;
  uint64_t n_reads;
};

// Helper for PM_LibraryStats that write the stats out in a nice table

void writeLibraryStats( ostream & out, const vec<PM_LibraryStats> & stats);


/******************************************************************************
 *
 *            PairsManager
 *
 ******************************************************************************/
class PairsManager {
  
private:
  
  // -------- MEMBER VARIABLES
  
  // Number of reads/pairs in the dataset.
  longlong _n_reads, _n_pairs;
  // IDs of the paired reads.  These two vecs are matched up so that
  // read _ID1[i] is paired to read _ID2[i].  The reads in a pair are not
  // guaranteed to appear in any particular order.  Unpaired reads will not be
  // listed here.
  vec<longlong> _ID1, _ID2;
  // Library IDs of the pairs.  These are indices into _libs.
  vec<unsigned char> _lib_IDs;
  // Library information.  The PM_Library struct is defined above.
  vec<PM_Library> _libs;

  // Cached data.  These vectors are filled when certain functions are called
  // and their contents are returned as references.  They will be cleared if
  // pairs are added or removed.
  mutable vec<longlong> _pairs_index, _partners_index;
  
  // allow each pair to have different seps and devs
  vec<int> seps, sds;
public:

  // -------- CONSTRUCTORS
  
  explicit PairsManager( const longlong & n_reads = 0 )
  { _n_reads = n_reads; _n_pairs = 0; }

  // Merge two PairsManagers from different datasets.
  PairsManager( const PairsManager & pm1, const PairsManager & pm2 );
  
  // Constructor that reads new PairManager format file.
  explicit PairsManager( const String & pairs_file )
  { Read( pairs_file ); }

  // Constructor that reads old-style pairto file.
  PairsManager( const String & pairto_file, const longlong n_reads )
  { ReadFromPairtoFile( pairto_file, n_reads ); }
  
  
  // -------- I / O 
  
  // Write a PairsManager to a file. Wrapper around writeBinary method.
  // By convention, these files usually have the extension '.pairs'.
  void Write( const String & pairs_file ) const;

  // Read a PairsManager from a file. Wrapper around readBinary method.
  void Read( const String & pairs_file );
  
  // Write a PairsManager to a stream, using BinaryWriter class in BinaryStream.h.
  // By convention, these files usually have the extension '.pairs'.
  size_t writeBinary( BinaryWriter& writer ) const;

  // Read a PairsManager from a stream, using the BinaryReader class in BinaryStream.h.
  // By convention, these files usually have the extension '.pairs'.
  void readBinary( BinaryReader& reader );

  static size_t externalSizeof() { return 0; }

  // Load a PairsManager from an old-style pairto file (e.g., reads.pairto.)
  // You must specify the number of reads because that information is not
  // contained in the file.
  void ReadFromPairtoFile( const String & pairto_file, const longlong n_reads );

  // Read only the library information from a PairsManager file.
  // Faster than loading the entire object.
  friend void ReadPairsManagerLibInfo( const String & pairs_file, longlong & nreads,
				       vec<String> & lib_names,
				       vec<int> & lib_sep, vec<int> & lib_sd );

  // Read only the library information from a PairsManager file.
  // Faster than loading the entire object.
  friend void ReadPairsManagerLibInfo( const String & pairs_file, longlong & nreads,
				       vec<PM_LibraryStats> & stats); 

private:

  // Read a PairsManager header from a stream, using the BinaryReader class in BinaryStream.h.
  // By convention, these files usually have the extension '.pairs'.
  void readHeader( BinaryReader& reader );

public:  
  
  // -------- FORMATTED PRINTING

  // Write nicely formatted library information to an output stream.
  void printLibraryStats( ostream & out, const String bases_filename = "" ) const;
  
  // Write nicely formatted information about each read to an output stream.
  void printPairInfo( ostream & out ) const;
  
  void print( ostream & out ) const { printLibraryStats( out ); out << "\n\n"; printPairInfo( out ); out << "\n\n"; }
  
  
  // -------- MODIFIERS

  void changeLibrarySepSd( const size_t & lib_ID, const int & sep, const int & sd ){
    _libs[ lib_ID ]._sep  = sep; _libs[ lib_ID ]._sd = sd;
  }

  void changeLibraryName( const size_t & lib_ID, const String & name ){
    _libs[ lib_ID ]._name  = name;
  }

  // These functions all call clearCache( ).
  
  void clear(); // use this to free up memory
  
  void SetIDs( const longlong & pair_ID, const longlong & ID1, const longlong & ID2 )
  { _ID1[pair_ID] = ID1; _ID2[pair_ID] = ID2; clearCache( ); }
  // Append a read pair to the list.  May require the creation of a new library;
  // if so, use lib_name (or default) as this library's name.
  void addPair( const longlong ID1, const longlong ID2,
		const int sep, const int sd,
		const String& lib_name = "",
		const bool grow_read_count = false );
  // Append a read pair to the list using existing library.
  void addPairToLib( const longlong ID1, const longlong ID2,
		     const size_t lib_ID,
		     const bool grow_read_count = false );

  // Append PairsManager to this PairsManager.
  void Append( const PairsManager & pm );
  
  void addLibrary( const int & sep, const int & sd, const String& lib_name = "");

  // Remove a set of read pairs.  DOES NOT remove the reads themselves.
  void removePairs( const vec<Bool> & to_remove );
  
  // Remove a set of reads, along with their pairing information.
  // If allow_severance = false, then any read whose partner is being removed
  // MUST also be removed, or the code will crash.
  // If allow_severance = true, then such a read will simply become unpaired.
  // The flag has no effect on reads that are already unpaired.
  void removeReads( const vec<Bool> & reads_to_remove,
		    const bool & allow_severance = false );
  
  
  // -------- QUERY FUNCTIONS: DATASET SIZES
  size_t nReads( ) const { return _n_reads; }
  size_t nPairs( ) const { return _n_pairs; }
  size_t nLibraries( ) const { return _libs.size(); }
  bool empty( ) const { return (_n_pairs == 0); }
  
  
  // -------- QUERY FUNCTIONS: INDIVIDUAL PAIRS
  // Note that, for the sake of efficiency, we do not check index bounds
  // in these functions.
  
  pair<longlong,longlong> operator[]( const longlong & pair_ID ) const {
    return make_pair( _ID1[pair_ID], _ID2[pair_ID] );
  }
  longlong ID1( const longlong & pair_ID ) const { return _ID1[pair_ID]; }
  longlong ID2( const longlong & pair_ID ) const { return _ID2[pair_ID]; }
  
  // Get either ID1 or ID2 from this pair - whichever one is not given as input.
  longlong getPartner( const longlong & pair_ID, const longlong & read_ID ) const;
  
  int libraryID( const longlong & pair_ID ) const { return _lib_IDs[pair_ID]; }
  String libraryName( const longlong & pair_ID ) const { return _libs[ _lib_IDs[pair_ID] ]._name; }
  //int sep( const longlong & pair_ID ) const { return _libs[ _lib_IDs[pair_ID] ]._sep; }
  //int sd ( const longlong & pair_ID ) const { return _libs[ _lib_IDs[pair_ID] ]._sd; }


  // allow each pairs to have different sep and sd
  int sep( const longlong & pair_ID ) const { 
    if (seps.empty()) return _libs[ _lib_IDs[pair_ID] ]._sep;
    else return seps[pair_ID]; }
  int sd ( const longlong & pair_ID ) const { 
    if (sds.empty()) return _libs[ _lib_IDs[pair_ID] ]._sd; 
    else return  sds[pair_ID];}
  void AddSeps(const vec<int>& seps_){ seps = seps_; }
  void AddSDs(const vec<int>& sds_){ sds = sds_; }
  int sep2( const longlong & pair_ID ) const { return seps[pair_ID]; }
  int sd2 ( const longlong & pair_ID ) const { return sds[pair_ID]; }
  
  // Get an individual value from the pairs index or partner index.
  // Unpaired reads will return -1.
  longlong getPairID   ( const longlong & read_ID ) const;
  longlong getPartnerID( const longlong & read_ID ) const;
  bool isPaired  ( const longlong & read_ID ) const { return ( getPairID(read_ID) != -1 ); }
  bool isUnpaired( const longlong & read_ID ) const { return ( getPairID(read_ID) == -1 ); }
  

  void changeLibraryId( const longlong & pair_ID, const int & lib_ID ) 
  {
    _lib_IDs[pair_ID] = lib_ID;
  }
  
  
  // -------- QUERY FUNCTIONS: CACHED DATA STRUCTURES
  // Index lists are returned as const references, meaning that they point
  // to objects stored inside the PairsManager.  Unpaired reads are marked in
  // index lists with -1.
  // WARNING: These variables will be erased if the PairsManager goes out of
  // scope, or if a "modifier" function is called, clearing the cache.
  // WARNING: These functions are not thread-safe.  If you intend to use them
  // in a multi-threaded environment, you should call makeCache() in a single
  // thread first.
  
  // Return an index list of each read's pair ID.  The values in this vec can
  // be used as inputs to the above functions ID1, ID2, sep, sd.
  vec<longlong> const & getPairsIndex( ) const;
  // Return an index list of each read's partner.
  vec<longlong> const & getPartnersIndex( ) const;
  
  
  // -------- QUERY FUNCTIONS: LIBRARY INFO
  
  // Functions to get info about all libraries at once.
  vec<String> getLibraryNames( ) const;
  vec<size_t> getLibrarySizes( ) const; // "size" = number of pairs
  vec<PM_LibraryStats> getLibraryStats( const String bases_filename = "" ) const;

  // Functions to get info about a library, given its ID in this PairsManager.
  String getLibraryName( const int & i ) const { return _libs[i]._name; }
  int getLibrarySep( const int & i ) const { return _libs[i]._sep; }
  int getLibrarySD ( const int & i ) const { return _libs[i]._sd; }
  void setLibrarySep( const int & i, int sep ) { _libs[i]._sep = sep; }
  void setLibrarySD ( const int & i, int sd ) { _libs[i]._sd = sd; }

  // Functions to get a library ID (in this PairsManager) from library info.
  // These functions return -1 if no library fits the description.
  int libraryID( const int & sep, const int & sd ) const;
  int libraryID( const String & lib_name ) const;
  
  
  // -------- MISCELLANEOUS
  
  // Convert these pairs into a vec<read_pairing> for use in legacy code.
  vec<read_pairing> convert_to_read_pairings( ) const;
  
  // Equality operator.  Requires ALL data to be the same.
  friend bool operator==( const PairsManager & pm1, const PairsManager& pm2 );
  
private:

  // -------- CACHING FUNCTIONS
  // These are declared const because the only variables they modify are cache
  // variables, which are mutable.
  // WARNING: These functions are non-thread-safe.  If you use a function that
  // calls make...IndexCache (e.g., getPairID, getPartnersIndex, etc.) from
  // inside a multithreaded function, it may cause a seg fault.  In this
  // situation, you should call makeCache( ) in a single thread first.
  
  void makePairsIndexCache( ) const;
  void makePartnersIndexCache( ) const;
  void clearCache( ) const { Destroy( _pairs_index, _partners_index ); }
  
public:
  void makeCache( ) const { makePairsIndexCache( ); makePartnersIndexCache( ); }
  
};

SELF_SERIALIZABLE(PairsManager);

void ReadPairsManagerLibInfo( const String & pairs_file, longlong & nreads,
			      vec<String> & lib_names,
			      vec<int> & lib_sep, vec<int> & lib_sd );
void ReadPairsManagerLibInfo( const String & pairs_file, longlong & nreads,
			      vec<PM_LibraryStats> & stats); 


#endif
