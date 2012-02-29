// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

#include "paths/VerifyGapLengths.h"
#include "paths/PathEmbedding.h"

void VerifyGapLengths( const MergedKmerPath& merger,
		       const KmerPath& p1,
		       const KmerPath& p2,
		       const NegativeGapValidator* ngv ) {

  const int p1start = merger.left_end.first;
  const int p2start = merger.left_end.second;

  const int p1end = merger.path.NSegments() - 1 - merger.right_end.first;
  const int p2end = merger.path.NSegments() - 1 - merger.right_end.second;

  const int left  = max( p1start, p2start );
  const int right = min( p1end, p2end );

  vec<PathEmbedding> junk, witness1, witness2;

  // When one read has multiple embeddings in the merger,
  // we can easily get a false positive.
  bool false_positive_likely =
    ( FindPathEmbeddings(p1, merger.path, junk, p1start, p1end) && junk.size()>1 ) ||
    ( FindPathEmbeddings(p2, merger.path, junk, p2start, p2end) && junk.size()>1 );

  for( int seg = left; seg <= right; seg++ )
    if( merger.path.isGap(seg) ) {
      KmerPath test_path = merger.path;
      const int low = merger.path.Minimum(seg), high=merger.path.Maximum(seg);

      // Is the minimal claimed gap length possible?
      test_path.SetGap( seg, low, low );
      if( !FindPathEmbeddings( p1, test_path, junk, p1start, p1end ) ) {
	cout << "\n\nGAP LENGTH ERROR!!!\n"
	     << "In segment " << seg << ", gap length " << low 
	     << " inconsistent with p1!" << endl;
	PRINT(merger.path); PRINT(p1); PRINT(p2);
	cout << endl;
      }

      if( !FindPathEmbeddings( p2, test_path, junk, p2start, p2end ) ) {
	cout << "\n\nGAP LENGTH ERROR!!!\n"
	     << "In segment " << seg << ", gap length " << low 
	     << " inconsistent with p2!" << endl;
	PRINT(merger.path); PRINT(p1); PRINT(p2);
	cout << endl;
      }
      
      // Is the maximal claimed gap length possible?
      test_path.SetGap( seg, high, high );
      if( !FindPathEmbeddings( p1, test_path, junk, p1start, p1end ) ) {
	cout << "\n\nGAP LENGTH ERROR!!!\n"
	     << "In segment " << seg << ", gap length " << high
	     << " inconsistent with p1!" << endl;
	PRINT(merger.path); PRINT(p1); PRINT(p2);
	cout << endl;
      }

      if( !FindPathEmbeddings( p2, test_path, junk, p2start, p2end ) ) {
	cout << "\n\nGAP LENGTH ERROR!!!\n"
	     << "In segment " << seg << ", gap length " << high
	     << " inconsistent with p2!" << endl;
	PRINT(merger.path); PRINT(p1); PRINT(p2);
	cout << endl;
      }
      
      // Is minimal claimed gap length really minimal?
      // Only check this if low is nonzero.
      if( low > 0 && (ngv == NULL ||
	  ngv->PossibleGap( test_path.Stop(seg-1), 
			    test_path.Start(seg+1),
			    low-1 )) ) {
	test_path.SetGap( seg, low-1, low-1 );
	if( FindPathEmbeddings( p1, test_path, witness1, p1start, p1end ) &&
	    FindPathEmbeddings( p2, test_path, witness2, p2start, p2end ) ) {
	  if( false_positive_likely )
	    cout << "\n\nPossible gap length error, but false positive more likely\n";
	  else
	    cout << "\n\nGAP LENGTH ERROR!!! -- probably\n";
	  cout << "In segment " << seg << ", gap length " << low-1
	       << " would work also!" << endl;
	  PRINT(merger.path); PRINT(test_path); PRINT(p1); PRINT(p2);
	  cout << endl;
	}
      }

      // Is maximal claimed gap length really maximal?
      if( ngv == NULL || ngv->PossibleGap( test_path.Stop(seg-1), 
					   test_path.Start(seg+1),
					   low+1 ) ) {
	test_path.SetGap( seg, high+1, high+1 );
	if( FindPathEmbeddings( p1, test_path, witness1, p1start, p1end ) &&
	    FindPathEmbeddings( p2, test_path, witness2, p2start, p2end ) ) {
	  if( false_positive_likely )
	    cout << "\n\nPossible gap length error, but false positive more likely\n";
	  else
	    cout << "\n\nGAP LENGTH ERROR!!! -- probably\n";
	  cout << "In segment " << seg << ", gap length " << high+1
	       << " would work also!" << endl;
	  PRINT(merger.path); PRINT(test_path); PRINT(p1); PRINT(p2);
	  cout << endl;
	}
      }
    }
}
