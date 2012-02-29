// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

#ifndef ALIGN_AND_MERGE
#define ALIGN_AND_MERGE

#include "CoreTools.h"
#include "Vec.h"

#include "paths/KmerPath.h"
// Just a forward declaration -- the .cc file has the include
class NegativeGapValidator;

/// The path merging tools return a MergedKmerPath for each merger.
/// This is mostly a KmerPath, but also includes some integers,
/// referring to segments in the merged path:
///  .given -- the interval aligning to the match originally passed 
///            to MergePaths
///  .left_end.{first,second} -- the interval containing the left endpoint
///            of the original paths
///  .right_end.{first,second} -- the interval containing the right endpoint
///            of the original paths, measured from the right end, so
///            path.NSegments()-1-right_end.first = where the first path stopped

// Each internal merging tool will only set the subset of these which it knows.


struct MergedKmerPath {
  // Empty constructor so I can make a vec of them.
  MergedKmerPath( ) { }
  MergedKmerPath( const KmerPath& rp ) { path = rp; }

  KmerPath path;
  pair<int,int> left_end;
  pair<int,int> right_end;

  int given;

  int longest_perfect_match;

  void flip() { 
    swap(left_end.first, left_end.second);
    swap(right_end.first, right_end.second);
  }

  friend ostream & operator<<( ostream &out, const MergedKmerPath &mrp )
  {
    for(int i=0; i<mrp.path.NSegments(); i++) {
      if(i==mrp.left_end.first) out << "|";
      if(i==mrp.left_end.second) out << "<";
      if(i==mrp.given) out << "*";
      out << mrp.path.Segment(i);
      if(i==mrp.given) out << "*";
      if(i==mrp.path.NSegments()-1-mrp.right_end.first) out << "|";
      if(i==mrp.path.NSegments()-1-mrp.right_end.second) out << ">";
    }
    
    return out;
  }

  void ColorPrint( ostream &out )
  {
    for(int i=0; i<path.NSegments(); i++) {
      if(i == left_end.second )
	out << START_MAGENTA << "<" << END_ESCAPE;
      if(i == left_end.first )
	out << START_BLUE << "<" << END_ESCAPE;
      if(i==given) out << "*";
      out << path.Segment(i);
      if(i==given) out << "*";
      if(i == path.NSegments()-1 - right_end.first )
	out << START_BLUE << ">" << END_ESCAPE;
      if(i == path.NSegments()-1 - right_end.second )
	out << START_MAGENTA << ">" << END_ESCAPE;
    }
    out << " (longest match=" << longest_perfect_match << ")" << endl;
  }

};

void MergePaths( const KmerPath& p1, const KmerPath& p2, 
		 int ind1, int ind2, vec<MergedKmerPath>& ans,
		 int min_perfect_match=1,
		 const NegativeGapValidator* ngv = NULL,
		 bool DEBUG_GAP_SIZES=false );



#endif
