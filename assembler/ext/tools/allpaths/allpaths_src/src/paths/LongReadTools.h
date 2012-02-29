///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef LONG_READ_TOOLS_H
#define LONG_READ_TOOLS_H

#include "Alignment.h"
#include "Basevector.h"
#include "CoreTools.h"
#include "PackAlign.h"
#include "graph/Digraph.h"
#include "kmers/KmerRecord.h"
#include "paths/AssemblyEdit.h"

// A gap patcher is defined by left and right ints u1 and u2, and a BaseVec r,
// whose left end aligns to u1 starting at tpos1 and whose right end aligns to u2
// ending at tpos2, where the positions are in terms of the BaseVecs associated
// to u1 and u2.  The idea is that the ends of r would have been trimmed back to
// places where it matches u1 and u2 along an L-mer.
//
//                                 tpos1             tpos2
//                                 *                 *
//    -----------------u1---------------    ---------------u2----------------
//                                 -------r-----------

class GapPatcher 
{
public:
  int t1, t2;
  BaseVec r;
  int rid;
  int tpos1, tpos2;
  int best_len;
  int expected_gap;
  
  GapPatcher() {}
  GapPatcher(const int t1, 
              const int t2, 
              const BaseVec& r,
              const int rid, 
              const int tpos1, 
              const int tpos2)
    : t1(t1), t2(t2), r(r), rid(rid), tpos1(tpos1), tpos2(tpos2),
      best_len(0), expected_gap(0) {}
  GapPatcher(const int t1, 
              const int t2, 
              const BaseVec& r,
              const int rid, 
              const int tpos1, 
              const int tpos2,
              const int best_len,
              const int expected_gap)
    : t1(t1), t2(t2), r(r), rid(rid), tpos1(tpos1), tpos2(tpos2),
      best_len(best_len), expected_gap(expected_gap) {}
  
  friend Bool operator==( const GapPatcher& p1, const GapPatcher& p2 )
  {   
    return (p1.t1    == p2.t1    && 
            p1.t2    == p2.t2    &&  
            p1.r     == p2.r     && 
            p1.rid   == p2.rid   && 
            p1.tpos1 == p2.tpos1 &&
            p1.tpos2 == p2.tpos2);
  }
  
  friend bool operator<( const GapPatcher& p1, const GapPatcher& p2 )
  {
    if (p1.t1    < p2.t1   ) return true;
    if (p1.t1    > p2.t1   ) return false;
    if (p1.t2    < p2.t2   ) return true;
    if (p1.t2    > p2.t2   ) return false;
    if (p1.r     < p2.r    ) return true;
    if (p1.r     > p2.r    ) return false;
    if (p1.rid   < p2.rid  ) return true;
    if (p1.rid   > p2.rid  ) return false;
    if (p1.tpos1 < p2.tpos1) return true;
    if (p1.tpos1 > p2.tpos1) return false;
    return (p1.tpos2 < p2.tpos2);
  }


  // ---- SELF_SERIALIZABLE method
  size_t writeBinary(BinaryWriter& writer) const
  {
    size_t len = 0;
    len += writer.write(rid);
    len += writer.write(t1);
    len += writer.write(tpos1);
    len += writer.write(t2);
    len += writer.write(tpos2);
    len += writer.write(r);
    return len;
  }
  // ---- SELF_SERIALIZABLE method
  void readBinary(BinaryReader& reader)
  {
    reader.read(&rid);
    reader.read(&t1);
    reader.read(&tpos1);
    reader.read(&t2);
    reader.read(&tpos2);
    reader.read(&r);
  }

  // ---- SELF_SERIALIZABLE method
  static size_t externalSizeof() { return 0; }
};

SELF_SERIALIZABLE(GapPatcher);






class GapPatcher0 
{
public:
  BaseVec r;
  int upos1, upos2;
  
  GapPatcher0() { }
  GapPatcher0(const BaseVec& r, const int upos1, const int upos2)
    : r(r), upos1(upos1), upos2(upos2) {}
  GapPatcher0(const GapPatcher & p)
    : r(p.r), upos1(p.tpos1), upos2(p.tpos2) {}
  
  friend bool operator==(const GapPatcher0& p1, const GapPatcher0& p2)
  { 
    return (p1.r     == p2.r     && 
            p1.upos1 == p2.upos1 && 
            p1.upos2 == p2.upos2);
  }
  
  friend bool operator!=(const GapPatcher0& p1, const GapPatcher0& p2)
  { return !(p1 == p2); }
  
  friend ostream& operator<<(ostream& out, const GapPatcher0& p)
  { return out << p.upos1 << ", " << p.upos2 << ", " << p.r.ToString( ); }
  
  friend bool operator<(const GapPatcher0& p1, const GapPatcher0& p2)
  {    
    if (p1.r     < p2.r)     return true;
    if (p1.r     > p2.r)     return false;
    if (p1.upos1 < p2.upos1) return true;
    if (p1.upos1 > p2.upos1) return false;
    return (p1.upos2 < p2.upos2);
  }
};



typedef GapPatcher0 gap_patcher0;

class align_data 
{
     public:

     align_data( ) { }

     align_data( const align& a, const int gid, const int pos2, const int Pos2,
          const int errors, const Bool fw )
          : a(a), gid(gid), pos2(pos2), Pos2(Pos2), errors(errors), fw(fw) { }

     friend Bool operator<( const align_data& a1, const align_data& a2 )
     {    if ( a1.errors < a2.errors ) return True;
          if ( a1.errors > a2.errors ) return False;
          if ( a1.gid < a2.gid ) return True;
          if ( a1.gid > a2.gid ) return False;
          if ( a1.pos2 < a2.pos2 ) return True;
          if ( a1.pos2 > a2.pos2 ) return False;
          if ( a1.Pos2 < a2.Pos2 ) return True;
          if ( a1.Pos2 > a2.Pos2 ) return False;
          return a1.fw < a2.fw;    }

     friend Bool operator==( const align_data& a1, const align_data& a2 )
     {    return a1.gid == a2.gid && a1.pos2 == a2.pos2 && a1.Pos2 == a2.Pos2
               && a1.errors == a2.errors && a1.fw == a2.fw;    }

     align a;
     int gid;
     int pos2, Pos2;
     int errors;
     Bool fw;

};

int DistX(const GapPatcher0& p1, 
          const GapPatcher0& p2,
          const BaseVec& b1, 
          const BaseVec& b2,
          alignment * p_al,
          double * timer);



void DeleteX(GapPatcher0 * p_p, 
             int * p_sum, 
             vec<alignment> * p_als,
             const vec<GapPatcher0>& P, 
             const int L,
             const BaseVec& b1, 
             const BaseVec& b2,
             double * timer);


void InsertX(GapPatcher0 * p_p, 
             int * p_sum, 
             vec<alignment> * p_als, 
             const vec<GapPatcher0>& P, 
             const int L,
             const BaseVec& b1, 
             const BaseVec& b2,
             double * timer);

void SubstituteX(GapPatcher0 * p_p, 
                 int * p_sum, 
                 vec<alignment> * p_als, 
                 const vec<GapPatcher0>& P, 
                 const int L,
                 const BaseVec& b1, 
                 const BaseVec& b2,
                 double * timer);

GapPatcher0 patcher_optimal_original2(const vec<GapPatcher0> & p0s,
                                      const size_t i_best, 
                                      const int L,
                                      const BaseVec & bv1,
                                      const BaseVec & bv2,
                                      double * timer = 0);


GapPatcher0 patcher_optimal_original(const vec<GapPatcher0> & p0s,
                                     const size_t i_best, 
                                     const int L,
                                     const BaseVec & bv1,
                                     const BaseVec & bv2,
                                     double * timer = 0);




int KmerId( const BaseVec& b, const int L, const int p );

Bool PerfectMatch( const BaseVec& b1, const BaseVec& b2,
     const int p1, const int p2, const int len );

// Define local alignments of a read to a sequence.  These are triples
//        (sequence id, position on read, offset)
// where the offset is the position on the sequence minus the position on the read. 
// An alignment is seeded on a matching L-mer whose flanks match "well enough".

void GetLocalAligns( const BaseVec& r, const BaseVecVec& U, 
     const vec< vec< pair<int,int> > >& Ulocs, const int L, const int flank,
     const int max_errs, vec< triple<int,int,int> >& aligns, ostream& out, 
     const Bool SHOW_ALIGNS = False, const Bool require_extend = False );

void GetGlobalAligns( const BaseVec& r, const BaseVecVec& U, 
     const vec< triple<int,int,int> >& aligns, const int bandwidth_div,
     vec<align>& aligns_a, const double sub_frac, const double ins_frac,
     const double del_frac );

// MakeAlignmentGraph.  The argument alignsx has the form
// ( u, { (upos1,rpos1), ..., (uposn,rposn) } )
// where 
// * u is an index in U
// * upos is a position on U[u]
// * rpos is a position on the sequence (e.g. a read) being aligned to U[u].
// The edge objects in G are offsets.

void MakeAlignmentGraph( const vec< pair< int, vec< pair<int,int> > > >& alignsx,
     const vecbasevector& U, const int verbosity, ostream& rout, digraphE<int>& G );

void AlignToRef( const vecbasevector& B, const String& run_dir,
     const String& data_dir, ostream& rout, const String& module );

void ConvertPathsIntoBases( const int L, const basevector& target,
     const digraphE<int>& G, const vec< vec<int> >& epaths,
     const vec<int>& vpaths, vec< pair< int, vec< pair<int,int> > > >& alignsx, 
     const vecbasevector& U, vec<basevector>& bpaths, const int verbosity, 
     ostream& rout, const String& run_dir, const String& data_dir,
     const String& module, const Bool simple = False );

Bool ConvertPathIntoBasesValid( const int L, const digraphE<int>& G,
     const vec<int>& p, vec< pair< int, vec< pair<int,int> > > >& alignsx,
     const vecbasevector& U );

int SmithWatFreeSym( const basevector& b1, const basevector& b2, align& a,
     const Bool penalize_left_gap = False, const Bool penalize_right_gap = False,
     unsigned int mismatch_penalty = 2, unsigned int gap_penalty = 3 );

class placementx {
     public:
     placementx( ) { }
     placementx( const int g, const int pos, const Bool fw )
          : g(g), pos(pos), fw(fw) { }
     int g;
     int pos;
     Bool fw;
};

vec<placementx> FindGenomicPlacements( basevector b, const int L,
     const vecbasevector& genome, const vec< vec< pair<int,int> > >& Glocs );


void patchers_gap_collect(const BaseVecVec & bvs,
                          const vec<GapPatcher> & patchers, 
                          const unsigned np_min, 
                          const int sz_padding_min,
                          vec<vec<GapPatcher> > * p_patchers_gap);

// type == 0   choose best patcher based on shortest patcher size 
// type == 1   choose best patcher based on median patcher size 
// type == 2   choose best patcher based on median gap size 

void patchers_first_guesses(const BaseVecVec & T,
                            const vec<vec<GapPatcher> > & patchers_gap,
                            const unsigned type,
                            const unsigned sz_patcher_min, // only for type 0
                            vec<unsigned> * p_ips_best);

// ---- sort gaps according to their predicted performance
//      time scales as  np . len^2 
// 
void indexes_gaps_for_parallel(const vec<vec<GapPatcher> >& patchers_gap, 
                               const vec<unsigned> & ips_best,
                               vec<unsigned> * p_is_gaps,
                               const bool verbose = false);

void patcher_short_print(const int i_gap,
                         const GapPatcher0 & p_opt, 
                         const BaseVec & bv1, 
                         const BaseVec & bv2);

class flank_datum {

     public:
     
     flank_datum( ) { }
     flank_datum( const int lshift, const int rshift, const int u1, const int u2,
          const int pos1, const int pos2, const basevector& target )
          : lshift(lshift), rshift(rshift), u1(u1), u2(u2), pos1(pos1), pos2(pos2),
          target(target) { }

     int lshift, rshift;
     int u1, u2;
     int pos1, pos2;
     basevector target;
     
};

template<int K> void CleanPatch(

     // inputs:

     const int L, const vecbasevector& U, const vec< vec< pair<int,int> > >& Ulocs, 
     const vec< triple<kmer<K>,int,int> >& kmers_plus,
     const fastavector& LEFT, const fastavector& RIGHT,

     // inputs and outputs:

     assembly_edit& e,

     // logging:

     ostringstream& outx, const Bool DIRECT, const int verbosity, 
     const String& data_dir, const String& run_dir );

template< int K > void GetFlankData( const fastavector& left, 
     const fastavector& right, const basevector& patch, 
     const vec< triple<kmer<K>,int,int> >& kmers_plus, 
     vec<flank_datum>& F, const int verbosity, ostream& rout, 
     const String& run_dir, const String& data_dir );

void AlignToTarget( 
     /* inputs */
     const basevector& target, const int K, const int L, const vecbasevector& U,
     const vec< vec< pair<int,int> > >& Ulocs,
     /* outputs */
     vec< pair< int, vec< pair<int,int> > > >& alignsx, 
     vec<double>& mismatch_ratesx, int& START, int& STOP,
     /* logging */
     const int verbosity, ostream& rout 
          );

void PickBestPath( const vec<basevector>& bpaths, const basevector& target,
     basevector& bpath, const int verbosity, ostream& rout, const String& run_dir,
     const String& data_dir );

#endif
