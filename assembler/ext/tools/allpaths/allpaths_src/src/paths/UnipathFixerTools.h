///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef UNIPATH_FIXER_TOOLS_H
#define UNIPATH_FIXER_TOOLS_H

#include "Basevector.h"
#include "CoreTools.h"
#include "Intvector.h"
#include "Qualvector.h"
#include "feudal/BinaryStream.h"
#include "kmers/KmerRecord.h"
#include "math/HoInterval.h"
#include "paths/Ulink.h"

// A chunk of work for the PatcherCottage.
struct PCottageJoinData
{
    int sep;
    int dev;
    bvec L;
    bvec R;
    vecbvec reads;
    vecqvec quals;
    vec< std::pair<int,int> > pairs;

    size_t writeBinary( BinaryWriter& writer ) const
    { size_t result = writer.write(sep);
      result += writer.write(dev);
      result += writer.write(L);
      result += writer.write(R);
      result += writer.write(reads);
      result += writer.write(quals);
      result += writer.write(pairs);
      return result; }

    void readBinary( BinaryReader& reader )
    { reader.read(&sep);
      reader.read(&dev);
      reader.read(&L);
      reader.read(&R);
      reader.read(&reads);
      reader.read(&quals);
      reader.read(&pairs); }

    static size_t externalSizeof() { return 0; }
    static size_t const HEADER = 0x415441444E494F4A; // "JOINDATA"
};
SELF_SERIALIZABLE(PCottageJoinData);

// Message to PatcherCottage describing which join to work on.
struct PCottageWhichJoin
{
    // which item to work on
    size_t itemNumber;

    // All the PCottageJoinData is written to a binary stream.  This member is
    // the offset for the BinaryReader of that PCottageJoinData -- seek to this
    // offset, and then read a PCottageJoinData struct.
    size_t offset;

    // It might be nice to rearrange the code so that these could be made a part
    // of the JoinData.
    int u1;
    int u2;
};
TRIVIALLY_SERIALIZABLE(PCottageWhichJoin);

// Message from PatcherCottage describing how it did on the join that was
// assigned to it.
struct PCottageResults
{
    size_t itemNumber;
    vecbvec reads;
    vec<int> startStop;
    String report; // it would be nice to break out the "perfect" indication
                   // so that we didn't have to do any string parsing

    size_t writeBinary( BinaryWriter& writer ) const
    { size_t result = writer.write(itemNumber);
      result += writer.write(reads);
      result += writer.write(startStop);
      result += writer.write(report);
      return result; }

    void readBinary( BinaryReader& reader )
    { reader.read(&itemNumber);
      reader.read(&reads);
      reader.read(&startStop);
      reader.read(&report); }

    static size_t externalSizeof() { return 0; }
};
SELF_SERIALIZABLE(PCottageResults);


// A segalign represents the alignment of a read to a unibase u (if fw), or an
// alignment of a read to the reverse complement of a unibase u (if !fw).  The
// conversion is (True, rid, rpos, rc(u), upos) --> (False, rid, rpos, u, upos).

class segalign {

     public:

     Bool fw;          // is read forward on unibase?
     uint64_t rid;     // read identifier
     int rpos;         // start position on read
     uint64_t u;       // unibase id;
     int upos;         // start position on unibase

     segalign( ) { }

     segalign( const Bool fw, const uint64_t rid, const int rpos, const uint u, 
          const int upos ) : fw(fw), rid(rid), rpos(rpos), u(u), upos(upos) { }

     // For sorting by u first, then upos.

     friend Bool operator<( const segalign& s1, const segalign& s2 )
     {    if ( s1.u < s2.u ) return True;
          if ( s1.u > s2.u ) return False;
          if ( s1.upos < s2.upos ) return True;
          if ( s1.upos > s2.upos ) return False;
          if ( s1.rid < s2.rid ) return True;
          if ( s1.rid > s2.rid ) return False;
          if ( s1.fw < s2.fw ) return True;
          if ( s1.fw > s2.fw ) return False;
          if ( s1.rpos < s2.rpos ) return True;
          return False;    }

};

// An opair is an oriented pair.  The first read is always "forward".

class opair {

     public:

     opair( ) { }
     opair( const int64_t id1, const int64_t id2, const int sep, const int dev )
          : id1(id1), id2(id2), sep(sep), dev(dev) { }

     int64_t id1, id2;
     int sep, dev;

     friend Bool operator==( const opair& x, const opair& y )
     {    return x.id1 == y.id1 && x.id2 == y.id2 && x.sep == y.sep
               && x.dev == y.dev;    }

     friend Bool operator<( const opair& x, const opair& y )
     {    if ( x.id1 < y.id1 ) return True;
          if ( x.id1 > y.id1 ) return False;
          if ( x.id2 < y.id2 ) return True;
          if ( x.id2 > y.id2 ) return False;
          if ( x.sep < y.sep ) return True;
          if ( x.sep > y.sep ) return False;
          if ( x.dev < y.dev ) return True;
          return False;    }

};

TRIVIALLY_SERIALIZABLE(opair);

class xalign {

     public:

     xalign( ) { }
     xalign( const int id1, const int id2, const int pos1, const int pos2 )
          : id1(id1), id2(id2), pos1(pos1), pos2(pos2) { }

     int32_t id1, id2;
     int32_t pos1, pos2;
     static vecbasevector* reads;

     friend Bool operator==( const xalign& a1, const xalign& a2 )
     {    return a1.id1 == a2.id1 && a1.id2 == a2.id2 
               && a1.pos1 == a2.pos1 && a1.pos2 == a2.pos2;    }

     friend Bool operator<( const xalign& a1, const xalign& a2 )
     {    if ( a1.id1 < a2.id1 ) return True;
          if ( a1.id1 > a2.id1 ) return False;
          if ( a1.id2 < a2.id2 ) return True;
          if ( a1.id2 > a2.id2 ) return False;
          if ( a1.pos1 < a2.pos1 ) return True;
          if ( a1.pos1 > a2.pos1 ) return False;
          if ( a1.pos2 < a2.pos2 ) return True;
          return False;    }

     int Len1( ) const { return (*reads)[id1].size( ); }
     int Len2( ) const { return (*reads)[id2].size( ); }

     int Offset( ) const { return pos1 - pos2; }

     // Return the overlap of the two reads in an alignment.

     int Overlap( ) const
     {    return IntervalOverlap( 0, Len1( ), Offset( ), Offset( ) + Len2( ) );    }

     int OverlapPlus( ) const
     {    return Overlap( ) + pos2;    }

     // RightExt2: return the amount that id2 extends id1 to the right, or zero if 
     // id1 is actually out front.  LeftExt1: the amount that id1 extends id2 to the
     // the left (or zero).

     int RightExt2( ) const { return Max( 0, Offset( ) + Len2( ) - Len1( ) ); }
     int LeftExt1( ) const { return Max( 0, Offset( ) ); }
     int LeftExt2( ) const { return Max( 0, -Offset( ) ); }

     // Return the interval of bases on id1 that are covered by the overlap with id2.

     ho_interval Cov1( ) const
     {    return ho_interval( pos1, pos1 + Overlap( ) );    }

     ho_interval Cov2( ) const
     {    return ho_interval( pos2, pos2 + Overlap( ) );    }

     int Errs( ) const
     {    int errs = 0;
          int start = pos1;
          int stop = Min( Len1( ), start + Len2( ) - pos2 );
          for ( int j = start; j < stop; j++ )
               if ( (*reads)[id1][j] != (*reads)[id2][j-pos1+pos2] ) ++errs;
          return errs;    }

     int Matches( ) const
     {    int matches = 0;
          int start = pos1;
          int stop = Min( Len1( ), start + Len2( ) - pos2 );
          for ( int j = start; j < stop; j++ )
               if ( (*reads)[id1][j] == (*reads)[id2][j-pos1+pos2] ) ++matches;
          return matches;    }

};

// An align_chain consists of a sequence of reads, together with a sequence of 
// alignments from each to the next, such that each read extends the previous one
// to the right by at least one base.

class align_chain {

     public:

     align_chain( ) { }
     align_chain( const vec<int>& ids, const vec<int>& pos1, const vec<int>& pos2,
          const int len )
          : ids(ids), pos1(pos1), pos2(pos2), len(len) { }
     align_chain( const int id ) 
     {    ids.push_back(id); 
          len = (*reads)[id].size( );    }

     static vecbasevector* reads;
     vec<int> ids;
     vec<int> pos1, pos2;

     int Len( ) const { return len; }

     void AddAlign( const xalign& a )
     {    ids.push_back(a.id2), pos1.push_back(a.pos1), pos2.push_back(a.pos2);    
          len +=  a.RightExt2( );    }

     basevector Seq( ) const
     {    ForceAssert( ids.nonempty( ) );
          basevector s = (*reads)[ ids[0] ];
          for ( int j = 0; j < ids.isize( ) - 1; j++ )
          {    int id1 = ids[j], id2 = ids[j+1];
               xalign a( id1, id2, pos1[j], pos2[j] );
               int rext = a.RightExt2( );
               if ( rext > 0 )
               {    basevector b;
                    const basevector& R = (*reads)[id2];
                    b.SetToSubOf( R, R.isize( ) - rext, rext );
                    s = Cat( s, b );    }    }
          return s;    }

     private:

     int len;

};

void AlignReadsToUnipaths( const String& run_dir, const String& jump_reads,
     const String& frag_reads, const String& frag_reads_edit, const Bool USE_JUMPS,
     const int MAX_PLACEMENTS, const String& unifile,
     vec< triple<int64_t,int64_t,int> >& ALIGNS,
     vec< triple<int64_t,int64_t,int> >& JALIGNS, 
     const String& checkpoint_head = "" );

void ExtendAligns( const int K, const vec< triple<int64_t,int64_t,int> >& ALIGNS,
     const vecbasevector& reads, const vecqualvector& quals,
     const vecbasevector& unibases, const vec< vec<int> > & nexts,
     const Bool print_alignments, const Bool print_segments, VecIntVec* calls,
     VecIntVec* qcalls, int64_t& total, int64_t& qtotal, int64_t& qgood,
     vec<segalign>& SEGS, const Bool off_end_ok, const Bool get_calls = True );

void UnibaseCopyNumbersFromSegAligns( const int K, const vecbasevector& unibases,
     const vec<int>& to_rc, const int PLOIDY, const String run_dir, 
     const String& FRAG_ALIGNS, const vec<segalign>& SEGS, 
     const vec<size_t>& S_START, vec<int>& CN );

void MakeUlinks( const String& run_dir, const String& TMP, const String& FRAG_READS,
     const String& JUMP_READS, const int K, const int PLOIDY, const vec<int>& CN, 
     const vec<int>& to_rc, vec<segalign>& SEGS, vec<segalign>& JSEGS,
     const vecbasevector& unibases, const vec<int>& innie_sep, 
     const vec<int>& innie_dev, const vec<double>& innie_percent, 
     const double min_innie_percent, const Bool log_links_detailed, 
     const int min_kmers, const int max_devs, vec<ulink_with_uids>& ulinks, 
     const Bool CHECKPOINT );

void BuildJoinData( const vec<int>& dead_fw, const vec<int>& dead_rc,
     const vec< pair<int,int> >& to_process_sepdev, const vecbasevector& unibases, 
     vec<int>& to_rc, const String& run_dir, const String& PATCHDIR,
     const String& FRAG_READS, const String& JUMP_READS, vec<segalign>& SEGS,
     vec<segalign>& JSEGS, vec<size_t>& S_START,
     vec<size_t>& JS_START, const vec<int>& innie_sep, const vec<int>& innie_dev, 
     const vec<double>& innie_percent, const double min_innie_percent,
     vec<size_t>& join_data_offsets, const vec< pair<int,int> >& joiners,
     const String& JOINDATA_file );

void MakeJoins( const vec< pair<int,int> >& Xjoiners, const vecbasevector& X,
     const String& progname, const int K, const String& data_dir,
     const String& run_dir, const String& PATCHDIR, int& attempted_joins,
     vec<Bool>& joined, int& npatches, vec<Bool>& perfect,
     bvec3& new_stuff, vec< vec<int> >& start_stop,
     vec<String>& reports, unsigned num_procs, const vec<int>& joins,
     const String& LOG, const Bool log_attempted_joins, const Bool log_action, 
     const int MAX_PATHS, const int MAX_PATH_ITERATIONS, const int MAX_EDGES, 
     const int MAX_READS, const int MAX_JOINS, const int MIN_OVERLAP_END,
     const String& JOINDATA_file, const vec<size_t>& join_data_offsets,
     const vec< pair<int,int> >& LR_to_process );

// BuildAlignments.  Find gap-free alignments between reads.  Subject to seed, 
// min_align_length, and max_align_errors, we find them all.  The computational 
// performance of this code might be improved in a couple of ways:
// - Reduce sizes of integers to minimum needed.  Would need to templatize to do 
//   this.
// - Don't create xaligns at all.

void BuildAlignments( const vecbasevector& reads, const int seed, 
     const int min_align_length, const double max_align_errors, 
     vec<xalign>& xaligns );

Bool Consistent( const align_chain& c, const xalign& a, const vec<xalign>& xaligns,
     const vec<int>& xaligns_index );

void GetRights( const align_chain& c, const vec<xalign>& xaligns,
     const vec<int>& xaligns_index, const vec<Bool>& good_right,
     vec<align_chain>& c_next, const int nreads, const int MIN_OVERLAP_END );

void CombinePairs( vecbasevector& reads, vecqualvector& quals, 
     const vec< pair<int,int> >& pairs, const double score_prox, 
     const int min_score, const int min_qdiff, const Bool verbosity, 
     String& report );

// Error correct each read using the reads that are aligned to it.  We do this by
// forming the pile of reads that align to a given read, and scoring each of its
// consecutive columns.  At present polymorphism is incorrectly handled.

void CorrectErrors( const vecbasevector& reads, const vecqualvector& quals,
     const vec<xalign>& xaligns, const int min_align_length,
     vecbasevector& reads_new, vecqualvector& quals_new,
     const Bool log_correct, String& report );

// Given unibases and a genome, report the number of genomic kmers missing from
// the unibases, the number of non-genomic kmers in the unibases, and the N50 
// unibase size.

template<int K> void UnibaseSummaryStats( const vecbasevector& unibases, 
     const vecbasevector& genome, vec<basevector>& missing_genomic_kmers,
     int64_t& non_genomic_kmers, int64_t& N50_unibase )
{    vec<int> unibase_sizes;
     for ( size_t i = 0; i < unibases.size( ); i++ )
          unibase_sizes.push_back( unibases[i].size( ) );
     Sort(unibase_sizes);
     N50_unibase = N50(unibase_sizes);
     missing_genomic_kmers.clear( );
     non_genomic_kmers = 0;
     int64_t nuni = unibases.size( ), nbases = 0, start = 0;
     vecbasevector F(unibases);
     F.Append(genome);
     vec<int64_t> starts;
     for ( size_t id = 0; id < F.size( ); id++ )
     {    starts.push_back(start);
          if ( F[id].isize( ) >= K )
          {    start += F[id].isize( ) - K + 1;
               nbases += F[id].isize( ) - K + 1;    }    }
     vec< kmer_record<K,2> > recs(nbases);
     #pragma omp parallel for 
     for ( size_t id = 0; id < F.size( ); id++ )
     {    const basevector& R2 = F[id];
          kmer_record<K,2> rec;
          basevector bi(K), bir(K);
          for ( int j = 0; j <= R2.isize( ) - K; j++ )
          {    bi.SetToSubOf( R2, j, K );
               bir.ReverseComplement(bi); 
               if ( bi < bir ) rec.Set( bi, id, j );
               else rec.Set( bir, id, -j-1 );
               recs[ starts[id] + j ] = rec;    }    }
     sort( recs.begin( ), recs.end( ), kmer_record<K,2>::cmp_basevector );
     for ( int i = 0; i < recs.isize( ); i++ )
     {    basevector b1, b2;
          recs[i].GetBasevector(b1);
          int j;
          for ( j = i + 1; j < recs.isize( ); j++ )
          {    recs[j].GetBasevector(b2);
               if ( b2 != b1 ) break;    }
          int genome_count = 0, unibases_count = 0;
          for ( int k = i; k < j; k++ )
          {    if ( recs[k].GetId( ) < nuni ) ++unibases_count;
               else ++genome_count;    }
          // The following fails, but shouldn't:
          // ForceAssert( unibases_count == 0 || unibases_count == 1 );
          if ( unibases_count == 0 ) 
          {    for ( int r = 0; r < genome_count; r++ )
                    missing_genomic_kmers.push_back(b1);    }
          if ( genome_count == 0 ) non_genomic_kmers += unibases_count;
          i = j - 1;    }    }

#endif
