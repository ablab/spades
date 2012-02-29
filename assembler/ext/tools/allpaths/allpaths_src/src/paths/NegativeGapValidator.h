// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

// A class to deal with K-mer read paths containing gaps of size <K.
// Uses a KmerBaseBroker to answer questions -- either pass it a
// pointer to one, or pass it the data to instantiate its own.

#ifndef NEGATIVEGAPVALIDATOR
#define NEGATIVEGAPVALIDATOR


#include "paths/KmerBaseBroker.h"
#include "paths/KmerPath.h"

// just need a forward declaration here; the .cc has the full include:
struct MergedKmerPath;

class NegativeGapValidator {

public:
  // constructor-destructor:
  NegativeGapValidator( const KmerBaseBroker* kbb ) :
    mp_kbb       ( kbb ), 
    I_own_my_kbb ( false ),
    K            ( kbb->GetK() ),
    verbosity      ( 0 )   { }

  NegativeGapValidator( String run_dir, int k ) :
    mp_kbb       ( new KmerBaseBroker(run_dir, k) ),
    I_own_my_kbb ( true ),
    K            ( k ),
    verbosity    ( 0 )   { }

  ~NegativeGapValidator( ) {
    if( I_own_my_kbb ) delete mp_kbb;
  }

  int GetK() const { return K; }

  void SetVerbosity( int v ) { verbosity = v; }


  // Just say yes or no
  bool Validate(const KmerPath& path, int seg_min=1, int seg_max=-2) const;

  // Fix if possible; say no otherwise
  // Note that this only adjusts gap lengths, so KmerPathLocs pointing
  // at kmers in the path remain valid.
  bool MakeValid(KmerPath& path, int seg_min=1, int seg_max=-2) const;

  // Apply the KBB's KmersBetween to each constant negative gap
  // Gaps may disappear, so KmerPathLocs may become invalid.
  bool FillConstantNegativeGaps(const KmerPath& path, KmerPath& ans,
				MergedKmerPath* mkp = NULL,
				int seg_min=1, int seg_max=-2 ) const;

  // Run MakeValid and then FillConstantNegativeGaps
  // KmerPathLocs may become invalid.
  bool MakeValidAndFillCNGs(KmerPath& path, 
			    MergedKmerPath* mkp = NULL,
			    int seg_min=1, int seg_max=-2 ) const {
    if( MakeValid( path, seg_min, seg_max ) ) {
      KmerPath ans;
      FillConstantNegativeGaps( path, ans, mkp, seg_min, seg_max );
      path = ans;
      return true;
    }
    else
      return false;
  }

  // Questions we just pass along to the KBB.
  // Syntactic sugar: versions for gaps, not offsets
  bool PossibleOffset( longlong k1, longlong k2, int d ) const
  { return mp_kbb->PossibleOffset(k1, k2, d); }
  bool PossibleGap( longlong k1, longlong k2, int gap ) const
  { return mp_kbb->PossibleOffset(k1, k2, gap+1); }

  int MinOffset( longlong k1, longlong k2, int d_min ) const
  { return mp_kbb->MinOffset(k1, k2, d_min); }
  int MinGap( longlong k1, longlong k2, int gap_min ) const
  { return( mp_kbb->MinOffset(k1, k2, gap_min+1)-1 ); }

  int MaxOffset( longlong k1, longlong k2, int d_max ) const
  { return mp_kbb->MaxOffset(k1, k2, d_max); }
  int MaxGap( longlong k1, longlong k2, int gap_max ) const
  { return( mp_kbb->MaxOffset(k1, k2, gap_max+1)-1 ); }

  // Question: would it be legal to concatenate these two paths?
  // The paths may begin or end with a gap interval.
  bool LegalConcat( const KmerPath& p, const KmerPath& q ) const;


private:
  const KmerBaseBroker* mp_kbb;
  bool I_own_my_kbb;
  int K;
  int verbosity;
};



#endif
