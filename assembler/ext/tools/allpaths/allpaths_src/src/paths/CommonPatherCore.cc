///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "paths/KmerPath.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/CommonPatherCore.h"

static inline 
String Tag(String S = "CP") { return Date() + " (" + S + "): "; } 

void CommonPather(const int K, 
                  const vec<String>& readsIn, 
                  const vec<String>& pathsOut,
                  const int NUM_THREADS,
                  const String CHECKPOINT_HEAD,
                  Bool appendKSize, 
                  Bool verify)
{
  bool VERBOSE = true;
  ForceAssertEq(readsIn.isize(), pathsOut.isize());

  String KS = ToString( K );

  int readSets = readsIn.isize();
  vec<int> setEnd(readSets);
  vecbasevector all;

  // Load basevectors
  for (int i = 0; i < readSets; i++) {
    cout << Tag() << "Loading " << Basename(readsIn[i]) << endl;
    all.ReadAll( readsIn[i], True );
    setEnd[i] = all.size();
    ForceAssertLe((i == 0 ? 0 : setEnd[i - 1]), setEnd[i]);
  }

  if ( all.size( ) == 0 )
  {    cout << "Something has gone terribly wrong.  Zero sequences were passed "
            << "to this module, CommonPather." << endl;
       cout << "Sorry, abort." << endl;
       exit(1);    }

  // Build paths
  cout << Tag() << "Building paths... " << endl;
  vecKmerPath all_paths;
  ReadsToPathsCoreY(all, K, all_paths, 
                    readsIn[0] + ".CommonPather", 
                    NUM_THREADS, CHECKPOINT_HEAD, VERBOSE);

  // Write paths
  for (int i = 0; i < readSets; i++) {
    String filename = pathsOut[i] + (appendKSize ? ".k" + KS : "");
    cout << Tag() << "Writing " << Basename(filename) << endl;
    all_paths.WriteRange(filename, (i == 0 ? 0 : setEnd[i - 1]), setEnd[i] );
  }

  // Code to test the repathing was sucessful.
  if (verify) {
    cout << Tag() << "Verifying sanity of pathed files..." << endl;

    for (int i = 0; i < readSets; i++) {
      String pathFile = pathsOut[i] + (appendKSize ? ".k" + KS : "");
      ForceAssert( IsGoodFeudalFile(pathFile));
      ForceAssertEq( MastervecFileObjectCount( pathFile ),
		     MastervecFileObjectCount( readsIn[i] ) );
    }

    vecKmerPath newAll;
    for (int i = 0; i < readSets; i++) {
      String pathFile = pathsOut[i] + (appendKSize ? ".k" + KS : "");
      newAll.ReadAll( pathFile + KS, True );
    }
    ForceAssertEq( all_paths.size(), newAll.size() );
    longlong mismatches = 0;
    for ( size_t i = 0; i < all_paths.size(); i++ ) {
      KmerPath p1 = all_paths[i];
      KmerPath p2 = newAll[i];
      if ( p1 != p2 ) {
	cout << Tag() << "path " << i << " does not match!!" << endl;
	mismatches++;
      }
    }
    if ( mismatches == 0 )
      cout << Tag() << "verified readback of " << newAll.size() << " paths." << endl;
    else {
      PRINT( mismatches );
    }
  }

  cout << Date( ) << ": Done!" << endl;
}


void CommonPather(const int K, 
                  const String & dirIn,
                  const vec<const vecbvec*>& readsIn, 
                  const vec<vecKmerPath*>& pathsOut,
                  const int NUM_THREADS,
                  const String CHECKPOINT_HEAD)
{
  bool VERBOSE = true;
  ForceAssertEq(readsIn.isize(), pathsOut.isize());

  int readSets = readsIn.isize();
  vec<int> setEnd(readSets);

  // Compute size and reserve merged vecbasevector
  size_t size = 0;
  for (int i = 0; i < readSets; i++) {
    size += readsIn[i]->size();
  }
  vecbasevector all;
  all.reserve(size);

  // Merge vecbasevectors
  for (int i = 0; i < readSets; i++) {
    all.Append(*readsIn[i]);
    setEnd[i] = all.size();
    ForceAssertLe((i == 0 ? 0 : setEnd[i - 1]), setEnd[i]);
  }

  // Build paths
  cout << Tag() << "Building paths... " << endl;
  vecKmerPath all_paths;
  ReadsToPathsCoreY(all, K, all_paths, 
                    dirIn + "/CommonPather", 
                    NUM_THREADS, CHECKPOINT_HEAD, VERBOSE);

  // Split paths
  for (int i = 0; i < readSets; i++) {
    pathsOut[i]->clear();
    pathsOut[i]->Append(all_paths, (i == 0 ? 0 : setEnd[i - 1]), setEnd[i]);
  }

}



void CommonPather(const int K, 
                  const String& dirIn, 
                  const vec<String>& readsInHead,
                  const int NUM_THREADS,
                  const String CHECKPOINT_HEAD)
{
  int readSets = readsInHead.isize();
  vec<String> readsIn(readSets);
  vec<String> pathsOut(readSets);
  for (int i = 0; i < readSets; i++) {
    readsIn[i] = readsInHead[i] + ".fastb";
    pathsOut[i] = readsInHead[i] + ".paths";
  }
  CommonPather(K, dirIn,
               readsIn, pathsOut,
               NUM_THREADS, CHECKPOINT_HEAD);
}


void CommonPather(const int K, 
                  const String& dirIn, 
                  const vec<String>& readsIn,
                  const vec<String>& pathsOut, 
                  const int NUM_THREADS,
                  const String CHECKPOINT_HEAD)
{
  ForceAssertEq(readsIn.isize(), pathsOut.isize());

  int readSets = readsIn.isize();
  vec<String> readsInFull(readSets);
  vec<String> pathsOutFull(readSets);
  for (int i = 0; i < readSets; i++) {
    readsInFull[i] = dirIn + "/" + readsIn[i];
    pathsOutFull[i] = dirIn + "/" + pathsOut[i];
  }
  CommonPather(K,
               readsInFull, pathsOutFull, 
               NUM_THREADS, CHECKPOINT_HEAD);
}



void CommonPather(const int K, 
                  const String & dirIn,
                  const vec<vecbvec*>& readsIn, 
                  const vec<vecKmerPath*>& pathsOut,
                  const int NUM_THREADS,
                  const String CHECKPOINT_HEAD)
{
  // Copy and cast vec of non-const vecbvec pointers to const vecbvec pointers
  vec<const vecbvec*> constReadsIn(readsIn.begin(),readsIn.end());
  CommonPather(K, 
               dirIn, 
               constReadsIn, pathsOut, 
               NUM_THREADS, CHECKPOINT_HEAD);
}


void CommonPather(const int K,
                  const String & dirIn,
                  const vecbvec& readsIn1, 
                  const vecbvec& readsIn2, 
                  const vecbvec& readsIn3,
                  vecKmerPath& pathsOut1, 
                  vecKmerPath& pathsOut2, 
                  vecKmerPath& pathsOut3,
                  const int NUM_THREADS,
                  const String CHECKPOINT_HEAD)
{
  
  vec<const vecbvec*> readsIn;
  vec<vecKmerPath*> pathsOut;

  readsIn.push_back(&readsIn1, &readsIn2);
  pathsOut.push_back(&pathsOut1, &pathsOut2);
  if (readsIn3.empty()) {
    readsIn.push_back(&readsIn3);
    pathsOut.push_back(&pathsOut3);
  }

  CommonPather(K, 
               dirIn, 
               readsIn, pathsOut, 
               NUM_THREADS, CHECKPOINT_HEAD);
   
}

void CommonPather(const int K,
                  const String & dirIn,
                  const vecbvec& readsIn1, 
                  const vecbvec& readsIn2,
                  vecKmerPath& pathsOut1, 
                  vecKmerPath& pathsOut2,
                  const int NUM_THREADS,
                  const String CHECKPOINT_HEAD)
{
  vecbvec dummy1;
  vecKmerPath dummy2;
  CommonPather(K, dirIn, 
               readsIn1, readsIn2, dummy1, pathsOut1, pathsOut2, dummy2,
	       NUM_THREADS, CHECKPOINT_HEAD);
}
