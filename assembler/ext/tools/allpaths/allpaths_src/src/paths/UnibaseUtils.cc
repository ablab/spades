/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "VecUtilities.h"
#include "paths/UnibaseUtils.h"
#include "kmers/KmerRecord.h"
#include "kmers/KmerShape.h"

/**
   Function: UnibaseInvolution

   For each <unibase>, identify its reverse complement.
 */
template < nbases_t K > 
void UnibaseInvolution( const vecbasevector& _unibases, vec< int >& _toRc ) {

  vecbasevector unibases( _unibases );
  
  vec<size_t> origUnibaseId( _unibases.size(), vec<size_t>::IDENTITY );
  
  unibases.SortSync( origUnibaseId );
  
  vec<size_t> toRc( unibases.size(), vec<size_t>::IDENTITY );
  
  vecbasevector unibases_rc( unibases );
  for ( size_t u = 0; u < unibases.size(); u++ )
    unibases_rc[ u ].ReverseComplement();
  
  unibases_rc.SortSync( toRc );
  for ( size_t u = 0; u < unibases.size(); u++ )
    if ( unibases[ u ] != unibases_rc[ u ] ) {
      PRINT( u );
      ForceAssertEq( unibases[ u ].ToString(), unibases_rc[ u ].ToString() );
    }
  
  /*
    {
    // check the results
    basevector u_rc;
    for ( size_t u = 0; u < unibases.size(); u++ ) {
    ForceAssertEq( unibases[ u ].ToString(), unibases_rc[ u ].ToString() );
    ForceAssertLt( toRc[ u ], unibases.size() );
    ForceAssertEq( toRc.size(), unibases.size() );
    ForceAssertEq( toRc[ toRc[ u ] ], u );
    u_rc.ReverseComplement( unibases[ u ] );
    ForceAssert( unibases[ toRc[ u ] ] == u_rc );
    }
    }
  */
  
  _toRc.resize( unibases.size(), -1 );
  for ( size_t i = 0; i < unibases.size(); i++ )
    _toRc[ origUnibaseId[ i ] ] = origUnibaseId[ toRc[ i ] ];
  
  /*
  // check the results
  basevector u_rc;
  for ( size_t u = 0; u < _unibases.size(); u++ ) {
  ForceAssertLe( 0, _toRc[ u ] );
  ForceAssertLt( static_cast<size_t>(_toRc[ u ]), _unibases.size() );
  ForceAssertEq( static_cast<size_t>(_toRc[ _toRc[ u ] ]), u );
  u_rc.ReverseComplement( _unibases[ u ] );
  ForceAssert( _unibases[ _toRc[ u ] ] == u_rc );
  }
  */

}

/**
   Function: UnibaseInvolution
   
   For each <unibase>, identify its reverse complement.
*/
void UnibaseInvolution( const vecbasevector& unibases, vec< int >& toRc, nbases_t K ) {
#define CASE(K) UnibaseInvolution<K>( unibases, toRc )
  DISPATCH_ON_K( K, CASE );
}





/**
   Function: BuildUnibaseAdjacencyGraph

 */
template < nbases_t K > 
void BuildUnibaseAdjacencyGraph( const vecbasevector& _unibases, digraph& _AG ) {

  _AG.Clear();
  vec< vec<int> > nexts;
  GetNexts( K, _unibases, nexts );
  int nuni = _unibases.size();

  vec< vec<int> > from(nuni), to(nuni);
  
  for ( size_t i = 0; i < nexts.size(); i++ ){
    for ( size_t j = 0; j < nexts[i].size(); j++ ){
      from[i].push_back(nexts[i][j]);
      to[ nexts[i][j] ].push_back( i );
    }
  }
  
  for ( size_t i = 0; i < from.size(); i++){
    Sort( from[i] );
    Sort( to[i] );
  }
  
  _AG.Initialize( from, to );
  ForceAssertEq( nuni, _AG.N() );

}

/**
   Function: BuildUnibaseAdjacencyGraph
   
*/
void BuildUnibaseAdjacencyGraph( const vecbasevector& unibases, digraph& AG, nbases_t K ) {
#define CASE2(K) BuildUnibaseAdjacencyGraph<K>( unibases, AG )
  DISPATCH_ON_K( K, CASE2 );
}


