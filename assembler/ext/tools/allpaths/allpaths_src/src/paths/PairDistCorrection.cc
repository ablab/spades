///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
#include "paths/PairDistCorrection.h"


// perform pair distance correction using the distribution obtained from IntDist.
// Note that the invariant separations are used in IntDist, while current MakeScaffold
// codes requires end-to-end separation. Therefore, read lengths are needed for the 
// conversion. 
//
void PairDistCorrectionFromIntDistOld(const String reads_head,  const PairsManager& pairs, const vec<alignlet>& aligns, const vec<int>& index, 
    vec<int>&seps, vec<int>& sds, bool verbose) {
  // get the average read length
  vec<PM_LibraryStats> stats = pairs.getLibraryStats(reads_head + ".fastb");
  // the original pair separations. Converted to invarant separation
  seps.resize(pairs.nPairs(),-1); 
  sds.resize(pairs.nPairs(),-1); 
  for( size_t i=0; i < pairs.nPairs(); i++ ) seps[i] =  pairs.sep(i);
  for( size_t i=0; i < pairs.nPairs(); i++ ) sds[i] =  pairs.sd(i);
  // read the distributions
  String file =  reads_head + ".distribs" ; 
  ForceAssert(IsRegularFile(file));
  vec<IntDistribution> distribs;
  BinaryReader::readFile(file.c_str(),&distribs);
  // initialize the probability functions for each library
  // note that the read lengths has to be subtracted 
  if (verbose) cout << Date() << ":     Loading library dist " <<endl;
  vec<ProbFuncIntDist> pfids;
  for(size_t i=0; i<distribs.size(); i++)
  {
    int libSep = pairs.getLibrarySep(i);
    int libSD =  pairs.getLibrarySD(i);
    pfids.push( ProbFuncIntDist(distribs[i],  -2*stats[i].mean_len, true ) ); // the distribution is shifted to match the end-to-end reads separation
  }
  if (verbose) cout << Date() << ":     End loading library dist " <<endl;
  PairDistCorrectionOld(pfids, pairs, aligns, index, seps, sds, verbose);
}


// perform pair distance correction using the distribution obtained from IntDist.
// Note that the invariant separations are used in IntDist, while current MakeScaffold
// codes requires end-to-end separation. Therefore, read lengths are needed for the 
// conversion. 
//
void PairDistCorrectionFromIntDist(const String reads_head,  const PairsManager& pairs, const vec<alignlet>& aligns, const vec<int>& index, 
    vec<int>&seps, vec<int>& sds, bool verbose) {
  // get the average read length
  vec<PM_LibraryStats> stats = pairs.getLibraryStats(reads_head + ".fastb");
  // the original pair separations. Converted to invarant separation
  seps.resize(pairs.nPairs(),-1); 
  sds.resize(pairs.nPairs(),-1); 
  // set seps[paidId] to be invarant separation
  int c = 0;
  for( size_t i=0; i < pairs.nPairs(); i++ ) {
    seps[i] =  pairs.sep(i) +  2 * stats[ pairs.libraryID(i) ].mean_len ;
  }
  for( size_t i=0; i < pairs.nPairs(); i++ ) sds[i] =  pairs.sd(i);
  // read the distributions
  String file =  reads_head + ".distribs" ; //  for distribution
  ForceAssert(IsRegularFile(file));
  vec<IntDistribution> distribs;
  BinaryReader::readFile(file.c_str(),&distribs);
  // initialize the probability functions for each library
  // note that the read lengths has to be subtracted 
  if (verbose) cout << Date() << ":     Loading library dist " <<endl;
  vec<ProbFuncIntDist> pfids;
  for(size_t i=0; i<distribs.size(); i++)
  {
    int libSep = pairs.getLibrarySep(i);
    int libSD =  pairs.getLibrarySD(i);
    pfids.push( ProbFuncIntDist(distribs[i] ) ); 
  }
  if (verbose) cout << Date() << ":     End loading library dist " <<endl;
  PairDistCorrection(pfids, pairs, aligns, index, seps, sds, verbose);
  // convert to end-to-end sep
  for( size_t i=0; i < pairs.nPairs(); i++ ) {
    seps[i] -= 2 * stats[ pairs.libraryID(i) ].mean_len ;
  }
}

// perform pair distance correction using the given distribution 
void PairDistCorrectionFromFile(const vec<String>& pairs_dist_files, const PairsManager& pairs, const vec<alignlet>& aligns, const vec<int>& index, 
    vec<int>&seps, vec<int>& sds, bool verbose) {
  seps.resize(pairs.nPairs(),-1); 
  sds.resize(pairs.nPairs(),-1); 
  for( size_t i=0; i < pairs.nPairs(); i++ ) seps[i] =  pairs.sep(i);
  for( size_t i=0; i < pairs.nPairs(); i++ ) sds[i] =  pairs.sd(i);
  // initialize the probability functions for each library
  vec<ProbFuncDist> pfds;
  for(size_t i=0; i<pairs_dist_files.size(); i++)
    pfds.push( ProbFuncDist(pairs_dist_files[i]) );
  PairDistCorrection(pfds, pairs, aligns, index, seps, sds, verbose);
}

