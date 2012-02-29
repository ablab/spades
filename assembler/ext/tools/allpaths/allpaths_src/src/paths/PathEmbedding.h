// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

#ifndef PATH_EMBEDDING_H
#define PATH_EMBEDDING_H

#include "CoreTools.h"
#include "Vec.h"

#include "paths/KmerPath.h"

/// A PathEmbedding is a data structure which represents the way a
/// subpath is embedded within a larger path.  That is, it is intended
/// for use only when every kmer of sub_path appears in super_path.
///
/// They are functionally read-only objects (with a private constructor).
/// To create a PathEmbedding, you should call FindPathEmbeddings()
/// (which may find more than one, in which case you will need to 
/// figure out which one you want), or one of the specialised generators,
/// IdentityEmbedding() and SubpathEmbedding().
///
/// INVALIDATION SEMANTICS: A PathEmbedding maintains a pointer to the
/// sub_path and super_path from which it is generated.  If these pointers
/// are invalidated, then the PathEmbedding's member functions can no longer
/// be used.  (They may give garbage answers or cause a segmentation fault.)
///
/// SPECIAL EMBEDDINGS: 
/// 1. An Identity embedding knows it's the identity because its sub_ and 
///    super_ pointers are the same.  It uses this knowledge for good.
/// 2. A call to SubpathEmbedding() creates a brand new read for the
///    subpath in question, which goes in the provided storage spot.
///    The caller must hold on to this, per invalidation semantics.

class PathEmbedding {
public:
  const KmerPath* GetSubPath() const { return sub_path; }
  const KmerPath* GetSuperPath() const { return super_path; }
  bool IsIdentity() const { return( sub_path == super_path ); }
  // theorem: any embedding of a path into itself must be the identity.
  
  int SubStartSegment() const { return block_start_on_super.front(); }
  int SubStopSegment() const { return block_stop_on_super.back(); }

  // Fast replacement for eg PushLoc(sub_path->End())
  KmerPathLoc SubStartLoc() const {
    return KmerPathLoc( super_path, block_start_on_super.front(),
			sub_path->Start(0) - 
			super_path->Start(block_start_on_super.front()) );
  }
  KmerPathLoc SubStopLoc() const {
    return KmerPathLoc( super_path, block_stop_on_super.back(),
			sub_path->Stop(sub_path->NSegments()-1) - 
			super_path->Start(block_stop_on_super.back()) );
  }



  // Push a KmerPathLoc on the sub_path through the embedding to the super_path
  KmerPathLoc PushLoc(const KmerPathLoc& sub_loc) const;

  // If a super_loc is in the image, find its preimage in sub_loc.
  // If the preimage is a gap, the Index() tells you which gap (ignore Loc()).
  KmerPathLoc Preimage(const KmerPathLoc& super_loc) const;

  // Quick check of whether super_loc is in the image at all:
  bool InImage(const KmerPathLoc& super_loc) const;

  // Push just a segment number from sub_path to super_path
  // No sense forcing the user to make a KmerPathLoc out of it just to do this.
  // This is well-defined but non-invertible.
  int PushSegment( int sub_seg ) const;


  // I guess these are unavoidable.  But I'm making them check path equality!
  void SetSubPath( const KmerPath* new_sub_path ) {
    ForceAssert( *sub_path == *new_sub_path );
    sub_path = new_sub_path;
  }
  void SetSuperPath( const KmerPath* new_super_path ) {
    ForceAssert( *super_path == *new_super_path );
    super_path = new_super_path;
  }


  friend bool operator<( const PathEmbedding& lhs, const PathEmbedding& rhs ) {
    return( lhs.block_start_on_super < lhs.block_start_on_super ||
	    (lhs.block_start_on_super == lhs.block_start_on_super &&
	     (*lhs.sub_path < *rhs.sub_path ||
	      (*lhs.sub_path == *rhs.sub_path && 
	       (*lhs.super_path < *rhs.super_path)))) );
    // The other two vec<int> members are just for efficiency --
    // they are guaranteed to be the same, if these all are.
  }

  friend bool operator==( const PathEmbedding& lhs, const PathEmbedding& rhs ){
    return( *lhs.sub_path == *rhs.sub_path &&
	    *lhs.super_path == *rhs.super_path &&
	    lhs.block_start_on_super == rhs.block_start_on_super );
    // The other two vec<int> members are just for efficiency --
    // they are guaranteed to be the same, if these all are.
  }


public:
  // For debugging mostly:
  friend ostream& operator<<(ostream& out, const PathEmbedding& e) {
    for( uint i=0; i<=e.sub_gap.size(); i++ ) {
      out << "sub block " << i << " runs from super segment "
	  << e.block_start_on_super[i] << " to "
	  << e.block_stop_on_super[i] << endl;
    }
    return out;
  }

private:
  // CONSTRUCTION:

  // PathEmbeddings are functionally read-only -- anything creating a new 
  // one should be a friend.  So the constructor is private.
  PathEmbedding(const KmerPath* sub, 
		const KmerPath* super) :
    sub_path(sub), super_path(super) {}
public:
  // Deep copy constructor -- copies referred-to paths into given new homes,
  // thereby guaranteeing the copied PathEmbedding remains valid.
  // SetSubPath and SetSuperPath obviate this, except for minor inefficiency
  PathEmbedding( const PathEmbedding& emb, 
		 KmerPath& new_sub_location, 
		 KmerPath& new_super_location) {
    ForceAssert( &new_sub_location != &new_super_location
		 || emb.IsIdentity() );
    *this = emb;   // copy over all the vec<int> data fields
    new_sub_location = *emb.GetSubPath();     // make a copy of emb->sub_path
    sub_path = &new_sub_location;             // and point to it
    new_super_location = *emb.GetSuperPath(); // make a copy of emb->super_path
    super_path = &new_super_location;         // and point to it
  }

  




  // Here's where they come from.  The return value is just !ans.empty()
  // If start_of_sub_on_super is specified, it's the segment of super that
  // the 0th segment of sub should align to.
  friend 
  bool FindPathEmbeddings( const KmerPath& sub, const KmerPath& super,
			   vec<PathEmbedding>& ans, 
			   int start_of_sub_on_super = -1,
			   int end_of_sub_on_super = -1 );

  friend 
  bool EmbeddingExists( const KmerPath& sub, const KmerPath& super,
			int start_of = -1, int end_of = -1 ) {
    vec<PathEmbedding> will_be_ignored;
    return( FindPathEmbeddings( sub, super, will_be_ignored,
				start_of, end_of ) );
  }


  friend PathEmbedding IdentityEmbedding( const KmerPath& path );

  // The subpath is stored in subpath_storage_spot
  friend PathEmbedding SubpathEmbedding( const KmerPath& path,
					 const KmerPathLoc& begin,
					 const KmerPathLoc& end,
					 KmerPath* subpath_storage_spot);


  // Compose the mappings from A->B and and B->C to a mapping of A->C.
  friend PathEmbedding ComposeEmbeddings( const PathEmbedding& embedAB,
                                          const PathEmbedding& embedBC );


  // This is a helper function for FindPathEmbeddings.
  // I doubt anything else should ever call it.
  friend 
  void RecursiveFindEmbeddings( const KmerPath& sub, const KmerPath& super,
				int sub_seg, int super_seg, 
				PathEmbedding& embedding,
				vec<PathEmbedding>& ans );
private:
  // Private helper functions:
  longlong first_kmer_in_sub_block( unsigned int b ) const
    { return sub_path->Start( b==0 ? 0 : sub_gap[b-1]+1 ); }
  longlong last_kmer_in_sub_block( unsigned int b ) const
    { return sub_path->Stop( b==sub_gap.size() ? sub_path->NSegments()-1
			     : sub_gap[b]-1 ); }


private:
  // Hold on to pointers to the paths, so that methods can return
  // useful information about them.  But this is dangerous, since
  // the caller might destroy (or edit) what these point to.
  const KmerPath* sub_path;
  const KmerPath* super_path;

  // These hold the number of the segment in super_path where
  // the ith gap-free-block of sub_path begins or ends.
  vec<int> block_start_on_super;
  vec<int> block_stop_on_super;  // For efficiency; technically not necessary
  // Also for efficiency, we record the gap segments in sub (in order):
  vec<int> sub_gap;
};

PathEmbedding IdentityEmbedding( const KmerPath& path );

PathEmbedding SubpathEmbedding( const KmerPath& path,
                                         const KmerPathLoc& begin,
                                         const KmerPathLoc& end,
                                         KmerPath* subpath_storage_spot);

bool FindPathEmbeddings( const KmerPath& sub, const KmerPath& super,
			 vec<PathEmbedding>& ans, 
			 int start_of_sub_on_super,
			 int end_of_sub_on_super );

bool EmbeddingExists( const KmerPath& sub, const KmerPath& super,
                        int start_of, int end_of );

#endif
