/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "FastIfstream.h" // fast_ifstream
#include "String.h"
#include "TokenizeString.h" // Tokenize
#include "ReadPairing.h" // read_pairing
#include "PairsManager.h"
#include "Basevector.h"
#include "feudal/BinaryStream.h"


// Write nicely formatted library information to an output stream.
void writeLibraryStats( ostream & out, const vec<PM_LibraryStats> & stats)
{

  if (stats.empty()) {
    out << "No libraries found!" << endl;
    return;
  }

  // Make a chart header line.
  vec< vec<String> > rows;
  vec<String> row1, row2;

  row1.push_back( "lib_name", "lib_ID", "sep", "dev" );
  row2.push_back( "--------", "------", "---", "---" );
  if (stats[0].n_reads != 0) {
    row1.push_back( "n_reads" );
    row2.push_back( "-------" );
    if (stats[0].n_bases != 0) {
      row1.push_back( "n_bases", "min_len", "max_len", "mean_len" );
      row2.push_back( "-------", "-------", "-------", "--------" );
    }
  }
  rows.push_back( row1, row2 );
  
  // For each library, make a line of info.
  for ( size_t i = 0; i < stats.size(); i++ ) {
    vec<String> row;
    if (stats[i].name == "Unpaired" )
      row.push_back( stats[i].name, "-", "-", "-" );
    else
      row.push_back( stats[i].name,
		     ToString(i),
		     ToString( stats[i].sep ),
		     ToString( stats[i].sd ) );
    if (stats[0].n_reads != 0) {
      row.push_back( ToString( stats[i].n_reads ) );
      if (stats[0].n_bases != 0)
	row.push_back( ToString(stats[i].n_bases),
		       ToString(stats[i].min_len),
		       ToString(stats[i].max_len),
		       ToString(stats[i].mean_len) );
    }
    rows.push_back(row);
  }
  
  // Convert the lines into a chart.
  PrintTabular( out, rows, 2, "lrrr" );
  out << endl;
}




// Write a PairsManager to a file. Wrapper around writeBinary method.
// By convention, these files usually have the extension '.pairs'.
void PairsManager::Write( const String & pairs_file ) const {
    BinaryWriter writer(pairs_file.c_str());
    writer.write(*this);
    writer.close();
}


// Read a PairsManager from a file. Wrapper around readBinary method.
void PairsManager::Read( const String & pairs_file ) {
  if ( !IsRegularFile( pairs_file ) )
    FatalErr( "ERROR: PairsManager::Read failed - can't find file " << pairs_file );
  BinaryReader reader(pairs_file.c_str());
  reader.read(this);
}


// Write a PairsManager to a stream, using BinaryWriter class in BinaryStream.h.
// By convention, these files usually have the extension '.pairs'.
size_t PairsManager::writeBinary( BinaryWriter& writer ) const {
  // File format versioning - some future proofing
  const int version = 1;
  size_t len = writer.write(version);

  // Write dataset size.
  len += writer.write(_n_reads);

  // Write library data.  In order to be consistent with file format version 1,
  // we have to write the library stats separate from the library names.
  // If we ever move to a version 2, we'll simplify this to:
  // writer.write( _libs );
  vec< pair<int,int> > lib_stats;
  for ( size_t i = 0; i < nLibraries(); i++ )
    lib_stats.push_back( make_pair( _libs[i]._sep, _libs[i]._sd ) );
  len += writer.write( lib_stats );
  len += writer.write( getLibraryNames() );

  // Write pair data.
  len += writer.write(_ID1);
  len += writer.write(_ID2);
  return len+writer.write(_lib_IDs);
}


// Read a PairsManager from a stream, using the BinaryReader class in BinaryStream.h.
// By convention, these files usually have the extension '.pairs'.
void PairsManager::readBinary( BinaryReader& reader ) {
  // Read pair header
  readHeader(reader);

  // Read pair data.
  reader.read(&_ID1);
  reader.read(&_ID2);
  reader.read(&_lib_IDs);

  // Update members
  _n_pairs = _ID1.size();

  // Sanity checks
  ForceAssertLe(_n_pairs * 2, _n_reads);
  ForceAssertEq(_ID1.size(), _ID2.size());
  
  clearCache( );
}

// Read a PairsManager header from a stream, using the BinaryReader class in BinaryStream.h.
void PairsManager::readHeader( BinaryReader& reader ) {
  // File format versioning - some future proofing
  int version;
  reader.read(&version);
  ForceAssertEq(version, 1);

  // Read dataset sizes.
  reader.read(&_n_reads);

  // Read library data.  In order to be consistent with file format version 1,
  // we have to read the library stats separate from the library names.
  // If we ever move to a version 2, we'll simplify this to:
  // reader.read( _libs );
  vec< pair<int,int> > lib_stats;
  vec<String> lib_names;
  reader.read(&lib_stats);
  reader.read(&lib_names);
  size_t n_libs = lib_stats.size();
  for ( size_t i = 0; i < n_libs; i++ )
    _libs.push_back( PM_Library( this, lib_stats[i].first, lib_stats[i].second, lib_names[i] ) );

  // Sanity checks
  ForceAssertEq(lib_stats.size(), lib_names.size());
}

// Load a PairsManager from an old-style pairto file (e.g., reads.pairto.)
// You must specify the number of reads because that information is not
// contained in the file.
void PairsManager::ReadFromPairtoFile( const String & pairto_file, const longlong n_reads )
{
  _n_reads = n_reads;
  _n_pairs = 0; // this is brought up later by calls to addPair( )
  
  if ( !IsRegularFile( pairto_file ) )
    FatalErr( "PairsManager: Can't find file: " << pairto_file );
  
  fast_ifstream in( pairto_file );
  
  // Read the first line of the pairto file to determine the number of pairs.
  String line;
  getline( in, line );
  longlong n_pairs_expected = line.Int( );
  ForceAssertGe( _n_pairs, 0 );
  ForceAssertLe( _n_pairs, _n_reads / 2 );
  
  // Clear and reserve space for pairs
  _ID1.clear();
  _ID2.clear();
  _lib_IDs.clear();
  _ID1.reserve( _n_pairs );
  _ID2.reserve( _n_pairs );
  _lib_IDs.reserve( _n_pairs );
  
  // Read the remaining lines of the pairto file to get the data for each pair.
  // If the number of lines (after the first one) is unequal to n_pairs_expected
  // there will be a failed assert.
  vec<String> tokens;
  for ( longlong i = 0; i < n_pairs_expected; i++ ) {
    getline( in, line );
    
    // Each line in the pairto file has 6 tokens, each with a specific meaning.
    // Get these tokens and create the pair.
    Tokenize( line, tokens );
    ForceAssertEq( tokens.size( ), 6u );
    addPair( tokens[0].Int( ), tokens[2].Int( ), tokens[1].Int( ), tokens[3].Int( ) );
  }
  
  // Verify that the file has no more lines (after a final newline.)
  getline( in, line );
  ForceAssert( in.fail( ) );
  
  clearCache( );
}



// Equality operator.  Requires ALL data to be the same.
bool operator==( const PairsManager & pm1, const PairsManager& pm2 )
{
  // For the sake of speed, we check the dataset sizes first, then
  // the smaller data structures, then the larger data structures.
  if ( pm1._n_reads != pm1._n_reads ) return false;
  if ( pm1._n_pairs != pm1._n_pairs ) return false;
  
  if ( pm1.nLibraries() != pm2.nLibraries() ) return false;
  for ( size_t i = 0; i < pm1.nLibraries(); i++ )
    if ( pm1._libs[i] != pm2._libs[i] ) return false;
  
  if ( pm1._lib_IDs != pm2._lib_IDs ) return false;
  if ( pm1._ID1 != pm2._ID1 ) return false;
  if ( pm1._ID2 != pm2._ID2 ) return false;
  
  return true;
}

  

// Write nicely formatted library information to an output stream.
void
PairsManager::printLibraryStats( ostream & out, const String bases_filename ) const
{
  vec<PM_LibraryStats> stats = getLibraryStats(bases_filename);

  writeLibraryStats(out, stats);
}


void
PairsManager::printPairInfo( ostream & out ) const
{
  // Print one line for each pair.
  cout << "  PAIR ID    READ 1    READ 2   LIB" << endl;
  cout << "  -------    ------    ------   ---" << endl;
  char line[40];
  for ( int64_t i = 0; i < _n_pairs; i++ ) {
    sprintf( line, "%9ld %9ld %9ld   %3d", i, _ID1[i], _ID2[i], int(_lib_IDs[i]) );
    out << line << endl;
  }
}



// Merge two PairsManagers from different datasets.
PairsManager::PairsManager( const PairsManager & pm1, const PairsManager & pm2 )
{
  // First, load the first input PairsManager (pm1) into this PairsManager.
  _n_reads = pm1._n_reads;
  _n_pairs = pm1._n_pairs;
  _ID1 = pm1._ID1;
  _ID2 = pm1._ID2;
  _lib_IDs = pm1._lib_IDs;
  _libs = pm1._libs;
  
  // Now append the second input PairsManager.
  Append( pm2 );
}



// Append PairsManager to this PairsManager.
void PairsManager::Append( const PairsManager & pm )
{
  // Bring in the pairs.
  // We must add them one at a time in order to merge the libraries properly.
  longlong _n_reads_orig = _n_reads;
  _n_reads += pm._n_reads;
  for ( longlong i = 0; i < pm._n_pairs; i++ ) {
    
    // Add an offset to the read IDs from pm2.
    longlong ID1 = pm.ID1( i ) + _n_reads_orig;
    longlong ID2 = pm.ID2( i ) + _n_reads_orig;
    int sep = pm.sep( i );
    int sd  = pm.sd ( i );
    
    // For all libraries in pm but not already in this PairsManager,
    // supply the correct library name.
    String lib_name = pm.libraryName( i );
    addPair( ID1, ID2, sep, sd, lib_name );
  }
}



// Clear this PairsManager and free up its memory.
void PairsManager::clear()
{
  clearCache();
  Destroy( _ID1, _ID2, _lib_IDs, _libs );
  _n_reads = _n_pairs = 0;
}



// Append a read pair to the list.  May require the creation of a new library.
void PairsManager::addPair( const longlong ID1, const longlong ID2,
			    const int sep, const int sd,
			    const String & lib_name,
			    const bool grow_read_count)
{
  // Check the validity of the read IDs.  Note that we do NOT verify that these
  // IDs have not already been assigned to other pairs.
  ForceAssertGe( ID1, 0 );
  ForceAssertGe( ID2, 0 );
  if (grow_read_count) {
    _n_reads = Max(ID1 + 1, ID2 + 1, _n_reads);
  } else {
    ForceAssertLt( ID1, _n_reads );
    ForceAssertLt( ID2, _n_reads );
  }
  _ID1.push_back( ID1 );
  _ID2.push_back( ID2 );
  
  // From the pair's separation and standard deviation, deduce the library ID.
  int lib_ID = -1;
  if (lib_name != "") {
    lib_ID = libraryID( lib_name );
    if (lib_ID != -1) // existing library (verify sep and sd are same)
      if ( _libs[lib_ID]._sep !=  sep || _libs[lib_ID]._sd !=  sd ) 
	FatalErr("Inconsistent sep/sd values for this library: " + lib_name);
  } else {
    lib_ID = libraryID( sep, sd );
  }
  
  // Create a new library if necessary (i.e., if no read with this sep/stdev
  // has already been seen.)  Use a default-generated library name.
  if ( lib_ID == -1 ) {
    lib_ID = _libs.size();
    PM_Library new_lib( this, sep, sd, lib_name );
    _libs.push_back( new_lib );
    
    ForceAssertLe( _libs.size(), numeric_limits<unsigned char>::max() );
  }
  _lib_IDs.push_back( lib_ID );
  
  _n_pairs++;
  ForceAssertLe( _n_pairs * 2, _n_reads );
  
  clearCache( );
}

// Append a read pair to the list using existing library.
void PairsManager::addPairToLib( const longlong ID1, const longlong ID2,
		   const size_t lib_ID,
		   const bool grow_read_count ) {

  // Check the validity of the read IDs.  Note that we do NOT verify that these
  // IDs have not already been assigned to other pairs.
  ForceAssertGe( ID1, 0 );
  ForceAssertGe( ID2, 0 );
  if (grow_read_count) {
    _n_reads = Max(ID1 + 1, ID2 + 1, _n_reads);
  } else {
    ForceAssertLt( ID1, _n_reads );
    ForceAssertLt( ID2, _n_reads );
  }
  _ID1.push_back( ID1 );
  _ID2.push_back( ID2 );
  
  // Check validity of the library
  ForceAssertLt( lib_ID, _libs.size() );

  _lib_IDs.push_back( lib_ID );
  
  _n_pairs++;

  ForceAssertLe( _n_pairs * 2, _n_reads );
  
  clearCache( );
}

void PairsManager::addLibrary( const int & sep, const int & sd, const String& name)
{
  // Create a new library
  PM_Library new_lib( this, sep, sd, name );
  _libs.push_back( new_lib );
  
  ForceAssertLe( _libs.size(), numeric_limits<unsigned char>::max() );
}



// Remove a set of read pairs.
void PairsManager::removePairs( const vec<Bool> & to_remove )
{
  ForceAssertEq( (longlong) to_remove.size(), _n_pairs );
  
  EraseIf( _ID1,     to_remove );
  EraseIf( _ID2,     to_remove );
  EraseIf( _lib_IDs, to_remove );
  _n_pairs = (longlong) _ID1.size( );
  clearCache( );
}



// Remove a set of reads, along with their pairing information.
// If allow_severance = false, then any read whose partner is being removed
// MUST also be removed (or the code will crash.)
// If allow_severance = true, then such a read will simply become unpaired.
// The flag has no effect on reads that are already unpaired.
void PairsManager::removeReads( const vec<Bool> & reads_to_remove, const bool & allow_severance )
{
  ForceAssertEq( (longlong) reads_to_remove.size(), _n_reads );
  
  // Create a mapping of old -> new read IDs.
  // Reads to be removed are mapped onto -1.
  vec<int64_t> ID_map( _n_reads, -1 );
  int64_t read_ID = 0;
  for ( int64_t i = 0; i < _n_reads; i++ ) {
    if ( reads_to_remove[i] ) continue;
    ID_map[i] = read_ID++;
  }
  
  // Reduce the number of reads.
  _n_reads = read_ID;
  
  
  // Re-number the reads within their pairings.
  vec<Bool> pairs_to_remove( _n_pairs, False );
  for ( int64_t i = 0; i < _n_pairs; i++ ) {
    _ID1[i] = ID_map[ _ID1[i] ];
    _ID2[i] = ID_map[ _ID2[i] ];
    
    // If either of these reads are being removed, we must remove the pair.
    bool remove1 = ( _ID1[i] == -1 );
    bool remove2 = ( _ID2[i] == -1 );
    
    if ( remove1 || remove2 ) {
      pairs_to_remove[i] = True;
      
      // If !allow_severance, and if only one of the reads in this pair is being
      // removed, then we complain.
      if ( !allow_severance && (!remove1 || !remove2) )
	FatalErr( "PairsManager::removeReads called with allow_severance = false, but a pair was severed." );
    }
  }
  
  // Remove all the pairs whose reads have been removed.
  removePairs( pairs_to_remove );
  clearCache();
}



// Get either ID1 or ID2 from this pair - whichever one is not given as input.
longlong PairsManager::getPartner( const longlong & pair_ID, const longlong & read_ID ) const
{
  if ( read_ID == _ID1[pair_ID] ) return _ID2[pair_ID];
  ForceAssertEq( read_ID, _ID2[pair_ID] ); // the input has to be equal to ONE of the reads
  return _ID1[pair_ID];
}



// Return an index list of each read's pair ID.  The values in this vec can
// be used as inputs to the functions ID1, ID2, sep, sd.
vec<longlong> const & PairsManager::getPairsIndex( ) const
{
  if ( _pairs_index.empty( ) ) makePairsIndexCache( );
  return _pairs_index;
}



// Return an index list of each read's partner.
vec<longlong> const & PairsManager::getPartnersIndex( ) const
{
  if ( _partners_index.empty( ) ) makePartnersIndexCache( );
  return _partners_index;
}



longlong PairsManager::getPairID( const longlong & read_ID ) const
{
  if ( _pairs_index.empty( ) ) makePairsIndexCache( );
  return _pairs_index[read_ID];
}



longlong PairsManager::getPartnerID( const longlong & read_ID ) const
{
  if ( _partners_index.empty( ) ) makePartnersIndexCache( );
  return _partners_index[read_ID];
}


// Returns the ID of the library with this sep/sd, or -1 if not found.
// Runtime is O(nLibraries()) - not efficient if many libraries (which is rare.)
int
PairsManager::libraryID( const int & sep, const int & sd ) const
{
  for ( size_t i = 0; i < nLibraries(); i++ )
    if ( _libs[i]._sep == sep && _libs[i]._sd == sd )
      return i;
  return -1; // -1 means "not found"
}


int PairsManager::libraryID( const String & lib_name ) const
{
  for ( size_t i = 0; i < nLibraries(); i++ )
    if ( _libs[i]._name == lib_name ) return i;
  
  return -1; // -1 means "not found"
}






vec<String>
PairsManager::getLibraryNames( ) const
{
  vec<String> names;
  for ( size_t i = 0; i < nLibraries(); i++ )
    names.push_back( getLibraryName(i) );
  return names;
}


vec<size_t>
PairsManager::getLibrarySizes( ) const
{
  // The size of each library is the number of pairs in it.
  vec<size_t> count( nLibraries(), 0 );
  for ( longlong i = 0; i < _n_pairs; i++ )
    count[ _lib_IDs[i] ]++;
  
  return count;
}

vec<PM_LibraryStats> 
PairsManager::getLibraryStats( const String bases_filename ) const
{
  vec<PM_LibraryStats> stats( nLibraries() + 1 ); // store unpaired info in last entry

  // populate basic library stats
  for ( size_t i = 0; i < nLibraries(); i++ ) {
    stats[i].name = getLibraryName(i);
    stats[i].sep = getLibrarySep(i);
    stats[i].sd = getLibrarySD(i);
  }
  stats[stats.size() - 1].name = "Unpaired";

  if (bases_filename == "") {
    for ( longlong i = 0; i < _n_pairs; i++ )
      stats[_lib_IDs[i]].n_reads += 2;
      stats[stats.size() - 1].n_reads = _n_reads - (_n_pairs * 2);
  } else {

    // Use VirtualMasterVec to avoid loading in vecbasevector
    typedef VirtualMasterVec<basevector> VVecBVec;
    VVecBVec bases((bases_filename).c_str());

    // compute read stats for each library
    size_t read_no = 0;
    for (VVecBVec::const_iterator basesItr = bases.begin();
	 basesItr != bases.end(); basesItr++) {
      longlong pairId = getPairID(read_no);
      int lib_no = (pairId != -1 ? libraryID(pairId) : nLibraries()); // store unpaired info in last entry
      uint32_t read_len= basesItr->size( );
      if (stats[lib_no].n_reads == 0) {
	stats[lib_no].min_len = stats[lib_no].max_len = stats[lib_no].n_bases = read_len;
      } else {
	if (read_len < stats[lib_no].min_len) stats[lib_no].min_len = read_len;
	if (read_len > stats[lib_no].max_len) stats[lib_no].max_len = read_len;
	stats[lib_no].n_bases += read_len;
      }
      stats[lib_no].n_reads++;
      read_no++;
    }
 
    // Compute additional library stats
    for ( size_t i = 0; i < stats.size(); i++ ) {
      stats[i].mean_len = (stats[i].n_reads == 0 ? 0 : stats[i].n_bases / stats[i].n_reads);
    }
  }
  
  if  (_n_reads ==  _n_pairs * 2) // No unpaired reads
    stats.resize(stats.size() - 1);

  return stats;
}


  
// Convert these pairs into a vec<read_pairing> for use in legacy code.
vec<read_pairing> PairsManager::convert_to_read_pairings( ) const
{
  vec<read_pairing> pairings( _n_pairs );
  
  for ( int i = 0; i < _n_pairs; i++ ) {
    pairings[i].id1 = _ID1[i];
    pairings[i].id2 = _ID2[i];
    pairings[i].sep = _libs[ _lib_IDs[i] ]._sep;
    pairings[i].sd  = _libs[ _lib_IDs[i] ]._sd;
    // Give the pairings generic values for t and weight.
    pairings[i].t = other;
    pairings[i].weight = 1;
  }
  
  return pairings;
}



// Fill _pairs_index.
void PairsManager::makePairsIndexCache( ) const
{
  if ( !_pairs_index.empty( ) ) return; // cache is already created
  _pairs_index.resize( _n_reads, -1 );
  
  for ( longlong i = 0; i < _n_pairs; i++ ) {
    // Verify that no read appears in multiple pairings.  (Note that we do
    // not generally make this check!)
    ForceAssertEq( _pairs_index[ _ID1[i] ], -1 );
    ForceAssertEq( _pairs_index[ _ID2[i] ], -1 );
    _pairs_index[ _ID1[i] ] = i;
    _pairs_index[ _ID2[i] ] = i;
  }
}



// Fill _partners_index.
void PairsManager::makePartnersIndexCache( ) const
{
  if ( !_partners_index.empty( ) ) return; // cache is already created
  _partners_index.resize( _n_reads, -1 );
  
  for ( longlong i = 0; i < _n_pairs; i++ ) {
    // Verify that no read appears in multiple pairings.  (Note that we do
    // not generally make this check!)
    ForceAssertEq( _partners_index[ _ID1[i] ], -1 );
    ForceAssertEq( _partners_index[ _ID2[i] ], -1 );
    _partners_index[ _ID1[i] ] = _ID2[i];
    _partners_index[ _ID2[i] ] = _ID1[i];
  }
}








// PM_Library constructor.  Provides a default name if none is given.
PM_Library::PM_Library( const PairsManager * pm, const int sep, const int sd, const String & name )
  : _pm ( pm ),
    _sep( sep ),
    _sd ( sd )
{
  ForceAssert( pm != NULL );
  
  if ( name != "" ) _name = name;
  else {
    // Default name is based on library ID, sep, and sd.
    _name = "LIBRARY" + ToString( pm->nLibraries() ) + "_sep" + ToString( _sep ) + "_sd" + ToString( sd );
  }
}




// PM_Library equality and inequality operators.  Require ALL data to be equal.
bool operator==( const PM_Library & lib1, const PM_Library& lib2 ) {
  if ( lib1._name != lib2._name ) return false;
  if ( lib1._sep  != lib2._sep  ) return false;
  if ( lib1._sd   != lib2._sd   ) return false;
  return true;
}
  
bool operator!=( const PM_Library & lib1, const PM_Library& lib2 ) {
  return !( lib1 == lib2 );
}


// Read only the library information from a PairsManager file.
// Faster than loading the entire object.
void ReadPairsManagerLibInfo( const String & pairs_file, longlong & nreads,
			      vec<String> & lib_names,
			      vec<int> & lib_sep, vec<int> & lib_sd ) {
  if ( !IsRegularFile( pairs_file ) )
    FatalErr( "ERROR: GetPairsManagerLibraryStats failed - can't find file " << pairs_file );
  BinaryReader reader(pairs_file.c_str());

  PairsManager pairs;

  // Read pair header
  pairs.readHeader(reader);

  nreads = pairs.nReads();
  size_t n_libs = pairs.nLibraries();
  lib_names.resize(n_libs);
  lib_sep.resize(n_libs);
  lib_sd.resize(n_libs);
  for ( size_t i = 0; i < n_libs; ++i ) {
    lib_names[i] = pairs.getLibraryName(i);
    lib_sep[i] = pairs.getLibrarySep(i);
    lib_sd[i] = pairs.getLibrarySD(i);
  }
}

// Read only the library information from a PairsManager file.
// Faster than loading the entire object.
void ReadPairsManagerLibInfo( const String & pairs_file, longlong & nreads, vec<PM_LibraryStats> & stats) {
  if ( !IsRegularFile( pairs_file ) )
    FatalErr( "ERROR: ReadPairsManagerLibInfo failed - can't find file " << pairs_file );

  BinaryReader reader(pairs_file.c_str());
  PairsManager pairs;

  // Read pair header
  pairs.readHeader(reader);

  nreads = pairs.nReads();
  size_t n_libs = pairs.nLibraries();

  stats = pairs.getLibraryStats();
  // Remove unpaired row - we don't know if there are any
  stats.resize(stats.size() - 1);
}

