// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

#include "paths/PathEmbedding.h"

// A few canonical embeddings:

PathEmbedding IdentityEmbedding( const KmerPath& path ) {
  PathEmbedding emb( &path, &path );
  emb.block_start_on_super.push_back(0);

  for( int seg=1; seg < path.NSegments()-1; seg++ )
    if( path.isGap(seg) ) {
      emb.block_stop_on_super.push_back(seg-1);
      emb.block_start_on_super.push_back(seg+1);
      emb.sub_gap.push_back(seg);
    }

  emb.block_stop_on_super.push_back(path.NSegments()-1);

  return emb;
}

// Compose embeddings of A into B and B into C into an embedding of A
// into C.

PathEmbedding ComposeEmbeddings( const PathEmbedding& embedAB,
                                 const PathEmbedding& embedBC )
{
  PathEmbedding embedAC( embedAB.GetSubPath(), embedBC.GetSuperPath() );

  for ( unsigned int i = 0; i < embedAB.block_start_on_super.size(); ++i )
  {
    KmerPathLoc blockStartOnB( embedAB.GetSuperPath(), embedAB.block_start_on_super[i] );
    blockStartOnB.SetKmer( embedAB.first_kmer_in_sub_block( i ) );

    KmerPathLoc blockStartOnC = embedBC.PushLoc( blockStartOnB );

    embedAC.block_start_on_super.push_back( blockStartOnC.GetIndex() );
  }

  for ( unsigned int i = 0; i < embedAB.block_stop_on_super.size(); ++i )
  {
    KmerPathLoc blockStopOnB( embedAB.GetSuperPath(), embedAB.block_stop_on_super[i] );
    blockStopOnB.SetKmer( embedAB.last_kmer_in_sub_block( i ) );

    KmerPathLoc blockStopOnC = embedBC.PushLoc( blockStopOnB );
    embedAC.block_stop_on_super.push_back( blockStopOnC.GetIndex() );
  }

  embedAC.sub_gap = embedAB.sub_gap;

  return embedAC;
}


/// This extracts a subpath of a path, and returns an embedding
/// of that subpath back into its parent path.  The subpath is
/// stored in the KmerPath* passed as an argument; the caller is
/// responsible for keeping this around as long as it is needed
/// for this or other related embeddings.

PathEmbedding SubpathEmbedding( const KmerPath& path,
				const KmerPathLoc& begin,
				const KmerPathLoc& end,
				KmerPath* subpath_storage_spot ) {  

  KmerPath& subpath = *subpath_storage_spot;
  subpath.Clear();
  path.CopySubpath( begin, end, subpath );

  PathEmbedding emb( subpath_storage_spot, &path );

  const int offset = begin.GetIndex();

  emb.block_start_on_super.push_back(offset);

  for( int seg=1; seg < subpath.NSegments()-1; seg++ )
    if( subpath.isGap(seg) ) {
      emb.block_stop_on_super.push_back(offset + seg - 1);
      emb.block_start_on_super.push_back(offset + seg + 1);
      emb.sub_gap.push_back(seg);
    }

  emb.block_stop_on_super.push_back(subpath.NSegments()-1 + offset);

  return emb;
}


void RecursiveFindEmbeddings( const KmerPath& sub, const KmerPath& super,
			      int sub_seg, int super_seg, 
			      PathEmbedding& embedding,
			      vec<PathEmbedding>& ans );

/// Find the embeddings of a subpath inside a superpath.
/// It's only an embedding if every kmer of the sub appears in
/// the super, with locations compatible with the sub's gaps.
/// For instance, after a merger, each of the merged paths is
/// embedded in the result.
///
/// It is possible (due to stretchy gaps or duplicated regions) 
/// that there will be multiple embeddings; this pushes all of
/// them onto the vector ans.  The return value is just !ans.empty().

bool FindPathEmbeddings( const KmerPath& sub, const KmerPath& super,
			 vec<PathEmbedding>& ans, 
			 int start_of_sub_on_super, int end_of_sub_on_super ) {
  ans.clear();

  // Definition: the empty path has no embeddings into anything.
  if( sub.IsEmpty() ) return false;

  PathEmbedding embedding( &sub, &super );

  longlong sub_first_kmer = sub.Start(0);
  
  if( start_of_sub_on_super != -1 ) {
    int seg = start_of_sub_on_super;
    if( super.Segment(seg).Contains(sub_first_kmer) )
      RecursiveFindEmbeddings( sub, super, 0, seg, embedding, ans );
  }
  else {
    for( int seg=0; seg < super.NSegments(); seg++ )
      if( super.isSeq(seg) && super.Segment(seg).Contains(sub_first_kmer) )
	RecursiveFindEmbeddings( sub, super, 0, seg, embedding, ans );
  }

  // Remove ones with wrong end_of_sub_on_super:
  if( end_of_sub_on_super != -1 ) {
    int good_count = 0;
    for( vec<PathEmbedding>::iterator i=ans.begin(); i!=ans.end(); i++ )
      if( i->block_stop_on_super.back() == end_of_sub_on_super )
	ans[good_count++] = *i;
    ans.erase( ans.begin()+good_count, ans.end() );
  }

  return !ans.empty();
}



// The PathEmbedding& embedding is used as a stack to record the
// current progress; when we find a success, we push a copy of it
// onto ans, and those make up our return values.

void RecursiveFindEmbeddings( const KmerPath& sub, const KmerPath& super,
			      int sub_seg, int super_seg, 
			      PathEmbedding& embedding,
			      vec<PathEmbedding>& ans ) {
  pair<KmerPathLoc,KmerPathLoc> locs = CreateAlignedLocs( sub, super, 
							  sub_seg, super_seg );
  KmerPathLoc& sub_loc = locs.first;
  KmerPathLoc& super_loc = locs.second;

  if( ! ScanRightPerfectMatch( sub_loc, super_loc ) // if improper or
      || (super_loc.atEnd() && !sub_loc.atEnd())    // sub extends past super or
      || (super_loc.GapToRight()                    // unexpected gap in super
	  && !sub_loc.GapToRight()
	  && !sub_loc.atEnd()) )
    return;

  // update the embedding to reflect the block we just aligned
  embedding.block_start_on_super.push_back(super_seg);
  embedding.block_stop_on_super.push_back(super_loc.GetIndex());

  // If we've reached the end of the subread, we win!
  // (Also in the weird case where sub ends with a gap.)
  if( sub_loc.atEnd()  ||  sub.NSegments() == sub_loc.GetIndex()+2 ) { 
    ans.push_back( embedding );
    // Pop the state of the current embedding:
    embedding.block_start_on_super.pop_back();
    embedding.block_stop_on_super.pop_back();
    return;
  }
  // Otherwise, there must be a gap to our right on sub.
  // Search super for the kmer just past the gap, and recurse.
  int sub_new_seg = sub_loc.GetIndex()+2;
  int seg = super_loc.GetIndex();
  embedding.sub_gap.push_back( sub_new_seg-1 );
  longlong sub_seek_kmer = sub.Start( sub_new_seg );

  int min_gap = sub.Minimum(sub_new_seg-1);
  int max_gap = sub.Maximum(sub_new_seg-1);

  // Special case where we stay in the same segment of super:
  if( super_loc.GetSegment().Contains( sub_seek_kmer ) &&
      sub_seek_kmer - sub_loc.GetKmer() - 1 >= min_gap &&
      sub_seek_kmer - sub_loc.GetKmer() - 1 <= max_gap ) {
    RecursiveFindEmbeddings( sub, super, sub_new_seg, seg,
			     embedding, ans );
  }
  // Now search forward for any other places where the gap might end
  int min_kmers_skipped = super_loc.Stop()-super_loc.GetKmer();
  int max_kmers_skipped = super_loc.Stop()-super_loc.GetKmer();

  for( seg++; seg < super.NSegments(); 
         min_kmers_skipped += super.MinLength(seg),
	 max_kmers_skipped += super.MaxLength(seg), 
	 seg++ ) {
    if( super.isSeq(seg) 
	&& super.Segment(seg).Contains( sub_seek_kmer )
	&& sub_seek_kmer - super.Start(seg) + max_kmers_skipped >= min_gap
	&& sub_seek_kmer - super.Start(seg) + min_kmers_skipped <= max_gap ) {
      RecursiveFindEmbeddings( sub, super, sub_new_seg, seg,
			       embedding, ans );
    }
  }

  // Now pop the state of the current embedding:
  embedding.block_start_on_super.pop_back();
  embedding.block_stop_on_super.pop_back();
  embedding.sub_gap.pop_back();
}



/// Take an RPL on sub_path, return the corresponding RPL on super_path.
/// This always exists.
//
// A KmerPathLoc holds a pointer to its path, which is why the
// PathEmbedding class keeps pointers to its sub and super paths.
// (Or the user could pass in sub and super each time... also seems bad)
// I wish I had a better way to design this.

KmerPathLoc PathEmbedding::PushLoc(const KmerPathLoc& sub_loc) const {

  if( IsIdentity() ) return sub_loc;

  // We only deal with sub_locs pointing to sequence.
  // I'm not sure what behavior would be desired for pushing gaps.
  if( !sub_loc.isSeq() )
    return KmerPathLoc( super_path, -1, -1 );

  // Figure out which gap-free-block the sub_loc points at:
  unsigned int block=0;
  while( block < sub_gap.size() && sub_gap[block] < sub_loc.GetIndex() )
    block++;

  // How many segments into the block are we?
  int seg_in_block = sub_loc.GetIndex();
  if( block > 0 ) 
    seg_in_block -= sub_gap[block-1] + 1;

  // What segment in the super path does that correspond to?
  int super_seg = block_start_on_super[block] + seg_in_block;

  // construct the return object
  KmerPathLoc super_loc( super_path, super_seg );
  super_loc.SetKmer( sub_loc.GetKmer() );  // this needs a valid super_path

  return( super_loc );
}

// Push just a segment number, not a whole loc:
int PathEmbedding::PushSegment( int sub_seg ) const {
  if( IsIdentity() ) return sub_seg;

  // We only deal with sub_locs pointing to sequence.
  // I'm not sure what behavior would be desired for pushing gaps.
  if( !sub_path->isSeq(sub_seg) )
    return( -1 );

  // Figure out which gap-free-block the sub_loc points at:
  unsigned int block=0;
  while( block < sub_gap.size() && sub_gap[block] < sub_seg )
    block++;

  if( block > 0 ) sub_seg -= sub_gap[block-1] + 1;
  
  return( block_start_on_super[block] + sub_seg );
}


/// Take an RPL on super_path, return its preimage on sub_path.
///
/// If the preimage is a gap, the returned RPL will point to the
/// gap segment; the loc is undefined.  Setting the loc to indicate
/// distance into the gap would require specifying distance from the
/// left or right edge, and how to deal with stretchy gaps on super_path.
/// If you want to get that information, you can take the returned RPL,
/// increment or decrement it one interval, push it to super_path,
/// and measure the minimum or maximum distance, as you desire.
///
/// If the super_loc is not in the image at all, it returns an
/// RPL with index set to -1.

KmerPathLoc PathEmbedding::Preimage(const KmerPathLoc& super_loc) const {
  
  if( IsIdentity() ) return super_loc;

  const int super_index = super_loc.GetIndex();
  const longlong super_kmer = super_loc.GetKmer();

  // fail quickly if super_loc is before image(sub) begins:
  if( super_index < block_start_on_super[0] ||
      (super_index == block_start_on_super[0] &&
       super_kmer < first_kmer_in_sub_block(0)) ) {
    return( KmerPathLoc( sub_path, -1, -1 ) );
  }

  unsigned int block = lower_bound( block_stop_on_super.begin(),
				    block_stop_on_super.end(),
				    super_index )- block_stop_on_super.begin();
  for( ; block < block_stop_on_super.size(); block++ ) {
    // Is the preimage in the preceding gap?
    if( block > 0 ) {
      if( (super_index > block_stop_on_super[block-1] ||
	   (super_index == block_stop_on_super[block-1] &&
	    super_kmer > last_kmer_in_sub_block(block-1)))
	  &&
	  (super_index < block_start_on_super[block] ||
	   (super_index == block_start_on_super[block] &&
	    super_kmer < first_kmer_in_sub_block(block))) ) {
	return( KmerPathLoc( sub_path, sub_gap[block-1] ) ); // no loc in gaps
      }
    }
    // Is the preimage in the current block?
    if( (super_index > block_start_on_super[block] ||
	 (super_index == block_start_on_super[block] &&
	  super_kmer >= first_kmer_in_sub_block(block)))
	&&
	(super_index < block_stop_on_super[block] ||
	 (super_index == block_stop_on_super[block] &&
	  super_kmer <= last_kmer_in_sub_block(block))) ) {
      int index = ( super_index - block_start_on_super[block]
                    + (block == 0 ? 0 : sub_gap[block-1]+1) );
      KmerPathLoc sub_loc(sub_path,index);
      sub_loc.SetKmer( super_kmer );
      return( sub_loc );
    }
  }
  // If we get here, super_loc is to the right of image(sub):
  return( KmerPathLoc( sub_path, -1, -1 ) );
}


bool PathEmbedding::InImage(const KmerPathLoc& super_loc) const {

  if( IsIdentity() ) return true;

  const int super_index = super_loc.GetIndex();
  const longlong super_kmer = super_loc.GetKmer();
  const unsigned int last_block = sub_gap.size();

  return( (super_index > block_start_on_super[0] ||
	   (super_index == block_start_on_super[0] &&
	    super_kmer >= first_kmer_in_sub_block(0)))
	  &&
	  (super_index < block_stop_on_super[last_block] ||
	   (super_index == block_stop_on_super[last_block] &&
	    super_kmer <= last_kmer_in_sub_block(last_block))) );
}
