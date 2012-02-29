///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef KMER_ALIGN_SET_H
#define KMER_ALIGN_SET_H

#include "Basevector.h"
#include "CoreTools.h"

// KmerAlignSet.  It represents a set of alignments of a read to unibases, although
// in principle the 'read' could be any sequence and the 'unibases' could be any set
// of sequences.
// - Each entry in the outer vec represents an 'alignment'.
// - First entry in the outer pair is the unibase id to which a read is aligned.
// - The inner pairs are (rpos,upos), presenting position on read, position on 
//   unibase.

class KmerAlignSet {

     public:

     KmerAlignSet( ) { }
     KmerAlignSet( const vec< pair< int, vec< pair<int,int> > > >& x ) : x_(x) { }
     KmerAlignSet( const basevector& x, const int K, 
          const vec< vec< pair<int,int> > >& Ulocs );

     int NAligns( ) const { return x_.size( ); }

     int U( const int i ) const { return x_[i].first; }

     int Count( const int i ) const { return x_[i].second.size( ); }

     int Rpos( const int i, const int j ) const { return x_[i].second[j].first; }
     int Upos( const int i, const int j ) const { return x_[i].second[j].second; }
     pair<int,int> RposUpos( const int i, const int j )
     {    return make_pair( Rpos(i,j), Upos(i,j) );    }
     int Offset( const int i, const int j ) const { return Upos(i,j) - Rpos(i,j); }

     void AddAlign( const int u, const vec< pair<int,int> >& rpos_upos )
     {    x_.push( u, rpos_upos );    }

     const vec< pair< int, vec< pair<int,int> > > >& X( ) const { return x_; }
     const pair< int, vec< pair<int,int> > >& X( const int n ) const 
     {    return x_[n];    }
     vec< pair< int, vec< pair<int,int> > > >& XMutable( ) { return x_; }

     private:

     vec< pair< int, vec< pair<int,int> > > > x_;

};

void ClusterAlignsNew( const KmerAlignSet& in, KmerAlignSet& out, const Bool clean,
     const int min_spread );

void ClusterAlignsOld( const KmerAlignSet& in, KmerAlignSet& out );

void KillInferiorClusters( KmerAlignSet& x );

void KillInferiorClustersNew( KmerAlignSet& x, const vecbasevector& unibases,
     const double min_ratio_to_kill = 5.0 );

#endif
