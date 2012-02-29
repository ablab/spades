// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

#ifndef MUTMERGRAPH
#define MUTMERGRAPH

#include <math.h>

#include "Basevector.h"
#include "CoreTools.h"
#include "math/Functions.h"
#include "pairwise_aligners/MaxMutmerFromMer.h"

// ================================================================================
//
// A mutmer_graph is a vector of mutmer_graph_node's.  Each node has a four-byte
// count, a four-byte next_node (see below), and room for a fixed number 
// (BLOCKS_PER_NODE) of mutmer_read_id's.  These come in small and big versions,
// depending on the the value of a template parameter I (1 or 2).  If I = 1, then
// the mutmer_read_id occupies 8 bytes, and contains the following data:
// 
// * read id of child (27 bits)
// * reverse complement? (1 bit)
// * pos1 = position on parent read (10 bits)
// * pos2 = position on child read (10 bits)
// * length = length of overlap (10 bits)
// * error count (5 bits) [Note: 31 means >= 31.]
//
// If I = 2, you get the big version, which occupies 20 bytes:
//
// * read id of child (31 bits)
// * reverse complement? (1 bit)
// * pos1 = position on parent read (32 bits)
// * pos2 = position on child read (32 bits)
// * length = length of overlap (32 bits)
// * error count (32 bits)
//
// The read id of the the parent is the index in the mutmer_graph of the
// mutmer_graph_node.
//
// If the reverse complement bit is set, then the region of overlap agrees as a
// reverse complement, e.g. as in:
//
//       ........ACGTT...    read 1      
//          .....AACGT.....  read 2
//
// Initially, next_node is set to -1.  Whenever a node becomes full 
// (count = BLOCKS_PER_NODE), a new node is allocated, and next_node is set to its
// index in overflow_.
//
// Note: Ony byte is enough for the count (so long as BLOCKS_PER_NODE is small).  So
// we can make count smaller if we need room for something else.
//
// This file should work for little endian architectures, but will need
// rewriting for the big endian case.
//
// ================================================================================

#ifndef Little_Endian
     #error MutmerGraph.h is designed for little endian architectures.
#endif

template<int I = 1> class mutmer_read_id {

     public:

     void Set( int read_id, int pos1, int pos2, int len, int errors, int rc )
     {    Assert( I == 1 || I == 2 );
          Assert( read_id >= 0 );
          Assert( rc == 0 || rc == 1 );
          if ( I == 1 )
	  {    Assert( read_id < ( 1 << 27 ) );
	       id_ = read_id | ((errors & 7) << 27) | (rc << 30);
               Assert( pos1 >= 0 && pos1 < 1024 );
               Assert( pos2 >= 0 && pos2 < 1024 );
               Assert( len >= 0 && len < 1024 );
               Assert( errors >= 0 && errors < 32 );
               mu_[0] = 
                    pos1 | (pos2 << 10) | (len << 20) | ((errors & 24) << 27);    }
          else if ( I == 2 )
          {    id_ = ( unsigned(read_id) << 1 ) ^ rc;
               mu_[0] = pos1;
               mu_[1] = pos2;
               mu_[2] = len;
               mu_[3] = errors;    }    }

     int ReadId( ) const
     {    if ( I == 1 ) return id_ & 0777777777;
          else return id_ >> 1;    }

     int Rc( ) const
     {    if ( I == 1 ) return (id_ & (1 << 30)) != 0;
          else return id_ & 1;    }

     unsigned int ReadIdRc( ) const
     {    if ( I == 1 ) return id_ & 010777777777;
          else return id_;    }

     void Unpack(int& pos1, int& pos2, int& len, int& e) const
     {    if ( I == 1 )
          {    pos1 = mu_[0] & 01777;
               pos2 = (mu_[0] >> 10) & 01777;
               len = (mu_[0] >> 20) & 01777;
               e = (((mu_[0] >> 30) & 3) << 3) | ((id_ >> 27) & 7);    }
          else if ( I == 2 )
          {    pos1 = mu_[0];
               pos2 = mu_[1];
               len = mu_[2];
               e = mu_[3];    }    }

     friend Bool operator<(const mutmer_read_id<I>& m1, const mutmer_read_id<I>& m2)
     {    return m1.ReadIdRc( ) < m2.ReadIdRc( );    }

     friend ostream& operator<<(ostream& s, const mutmer_read_id<I>& m)
     {    s << "read id: " << m.ReadId( );
          if ( m.Rc( ) ) s << " (rc)";
          int pos1, pos2, len, e;
          m.Unpack( pos1, pos2, len, e );
          s << ", pos1 = " << pos1 << ", pos2 = " << pos2 
               << ", length = " << len << ", errors = " << e << "\n";
          return s;    }

     unsigned int Id( ) const
     {    return id_;    }

     // THE FOLLOWING ARE ONLY VALID IF I = 2!!!

     int pos1( ) const { return mu_[0]; }
     int pos2( ) const { return mu_[1]; }
     int len( ) const { return mu_[2]; }
            
     private:

     unsigned int id_;
     unsigned int mu_[I*I];

};

template<int I, int BLOCKS_PER_NODE> class mutmer_graph_node {

     public:

     mutmer_graph_node( ) 
     {    count = 0;
          next_node = -1;   }

     void reset( )
     {    count = 0;
          next_node = -1;    }

     int count;
     int next_node;
     mutmer_read_id<I> contents[BLOCKS_PER_NODE];

};

inline int FirstBeforeSecond( const int& read_id1, const int& read_id2 )
{    if ( read_id1 % 1000 < read_id2 % 1000 ) return 1;
     if ( read_id1 % 1000 > read_id2 % 1000 ) return 0;
     return read_id1 < read_id2;    }

template<int I, int BLOCKS_PER_NODE> class mutmer_graph {

     public:

     mutmer_graph( ) { }

     mutmer_graph(int N)
     {    node_.resize(N);
          overflow_count_ = 0;
          unsigned int initial_overflow_size = Max( N/10, 2 );
          ForceAssert( initial_overflow_size < (((unsigned int)(1))<<31) );
          overflow_.resize(initial_overflow_size);    }

     void clear_and_resize(int N)
     {    node_.resize(N);
          for ( int i = 0; i < N; i++ )
               node_[i].reset( );
          overflow_count_ = 0;
          unsigned int initial_overflow_size = Max( N/10, 2 );
          ForceAssert( initial_overflow_size < (((unsigned int)(1))<<31) );
          overflow_.resize(initial_overflow_size);
          for ( unsigned int i = 0; i < overflow_.size( ); i++ )
               overflow_[i].reset( );    }

     mutmer_graph_node<I, BLOCKS_PER_NODE>* Overflow(unsigned int i)
     {    Assert( i < overflow_count_ );
          return &(overflow_[i]);    }

     vec<int> Counts( )
     {    vec<int> stats( node_.size( ) );
          for ( unsigned int i = 0; i < stats.size( ); i++ )
          {    stats[i] = node_[i].count;
               mutmer_graph_node<I, BLOCKS_PER_NODE>* n = &(node_[i]);
               while ( n->next_node >= 0 )
               {    n = Overflow( n->next_node );
                    stats[i] += n->count;    }    }
          return stats;    }

     int MaxCount() {
       int max_count = -1;
       for ( unsigned int i = 0; i < node_.size( ); i++ ) {
	 mutmer_graph_node<I, BLOCKS_PER_NODE>* n = &(node_[i]);
	 int total = n->count; 
	 while ( n->next_node >= 0 ) {
	   n = Overflow( n->next_node );
	   total += n->count;
	 }
	 if(total > max_count) max_count = total;
       }
       return max_count;
     }

     int OverflowCount( )
     {    return overflow_count_;    }

     void AllocateOverflowNode( )
     {    ++overflow_count_;
          if ( overflow_count_ >= overflow_.size( ) )
          {    unsigned int new_overflow_size = 
                    overflow_.size( ) + overflow_.size( )/2;
               Assert( new_overflow_size < (((unsigned int)(1))<<31) );
               overflow_.resize(new_overflow_size);
               // cout << "resizing overflow_ to " << overflow_.size( ) << "\n";   
                    }   }

     int All( int n, vec< mutmer_read_id<I> >& v )
     {    int count = 0;
          for ( int i = 0; i < node_[n].count; i++ )
               v[count++] = node_[n].contents[i];
          int next = node_[n].next_node;
          while ( next != -1 )
          {    mutmer_graph_node<I, BLOCKS_PER_NODE>* m = Overflow(next);
               for ( int i = 0; i < m->count; i++ )
                    v[count++] = m->contents[i];
               next = m->next_node;    }
          return count;    }

     void MergeKmer( int pos1, int pos2, int K, const int& read_id1, 
          const int& read_id2, const basevector& rd1, const basevector& rd2,
          Bool strict = False )
     {    
          Assert( pos1 != 0 && pos2 != 0 );

          // Set rc if exactly one of pos1, pos2 are negative.  Then 
          // set pos_i to abs(pos_i) - 1, for i = 1, 2.

          int rc = (pos1 < 0) ^ (pos2 < 0);
          if ( pos1 < 0 ) pos1 = -pos1;
          if ( pos2 < 0 ) pos2 = -pos2;
          --pos1;
          --pos2;

          Assert( pos1 >= 0 && pos2 >= 0 );

          int xpos1, xpos2, xk, xerrors;
          unsigned int read_id1_rc, read_id2_rc;
          if ( I == 2 ) 
          {    read_id1_rc = (read_id1 << 1) ^ rc;
               read_id2_rc = (read_id2 << 1) ^ rc;    }
          mutmer_graph_node<I, BLOCKS_PER_NODE>* n;
          if ( FirstBeforeSecond( read_id1, read_id2 ) ) 
          {    n = &(node_[read_id1]);
               if (rc) pos2 = rd2.size( ) - K - pos2;
               top1:
               for ( int i = 0; i < n->count; i++ )
               {    mutmer_read_id<I>& m = n->contents[i];
                    if ( I == 1 )
                    {    if ( m.ReadId( ) == read_id2 && rc == m.Rc( ) )
                         {    m.Unpack(xpos1, xpos2, xk, xerrors);
                              if ( xpos1 - xpos2 != pos1 - pos2 ) continue;
                              if ( pos1 >= xpos1 && pos1 + K <= xpos1 + xk ) 
                                   return;    }    }
                    else
                    {    if ( m.ReadIdRc( ) == read_id2_rc )
                         {    if ( m.pos1( ) - m.pos2( ) != pos1 - pos2 ) continue;
                              if ( pos1 >= m.pos1( ) && pos1 + K 
                                   <= m.pos1( ) + m.len( ) ) return;    }    }    }
               if ( n->next_node >= 0 )
               {    n = Overflow( n->next_node );
                    goto top1;    }
               int errors = 0;
               if (!rc) MaxMutmerFromMer( pos1, pos2, K, errors, rd1, rd2, strict );
               else MaxMutmerFromMerRev( pos1, pos2, K, errors, rd1, rd2, strict );
               mutmer_read_id<I>& m = n->contents[n->count];
               m.Set( read_id2, pos1, pos2, K, errors, rc );    }
          else
          {    n = &(node_[read_id2]);
               if (rc) pos1 = rd1.size( ) - K - pos1;
               top2:
               for ( int i = 0; i < n->count; i++ )
               {    mutmer_read_id<I>& m = n->contents[i];
                    if ( I == 1 )
                    {    if ( m.ReadId( ) == read_id1 && rc == m.Rc( ) )
                         {    m.Unpack(xpos1, xpos2, xk, xerrors);
                              if ( xpos1 - xpos2 != pos2 - pos1 ) continue;
                              if ( pos2 >= xpos1 && pos2 + K <= xpos1 + xk ) 
                                   return;    }    }
                    else
                    {    if ( m.ReadIdRc( ) == read_id1_rc )
                         {    if ( m.pos1( ) - m.pos2( ) != pos2 - pos1 ) continue;
                              if ( pos2 >= m.pos1( ) && pos2 + K 
                                   <= m.pos1( ) + m.len( ) ) return;    }    }    }
               if ( n->next_node >= 0 )
               {    n = Overflow( n->next_node );
                    goto top2;    }
               int errors = 0;
               if (!rc) MaxMutmerFromMer( pos2, pos1, K, errors, rd2, rd1, strict );
               else MaxMutmerFromMerRev( pos2, pos1, K, errors, rd2, rd1, strict );
               mutmer_read_id<I>& m = n->contents[n->count];
               m.Set( read_id1, pos2, pos1, K, errors, rc );    }
          ++n->count;
          if ( n->count == BLOCKS_PER_NODE )
          {    n->next_node = OverflowCount( );
               AllocateOverflowNode( );    }    }

     private:

     vec< mutmer_graph_node<I, BLOCKS_PER_NODE> > node_;
     unsigned int overflow_count_;
     vec< mutmer_graph_node<I, BLOCKS_PER_NODE> > overflow_;
};

#endif
