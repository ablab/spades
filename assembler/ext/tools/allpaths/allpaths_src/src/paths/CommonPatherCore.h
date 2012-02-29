/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/// CommonPather
/// 
/// Generalized read pathing tool that paths multiple sets of reads into the
/// same kmer space. Given a list of fastb files, it will produce a corresponding
/// set of paths files which share the same kmer space. Alternatively, it will work
/// directly on a set of vecbasevector, returning the corresponding vecKmerPaths.

#ifndef COMMON_PATHER_CORE_H
#define COMMON_PATHER_CORE_H

#include "CoreTools.h"
#include "Basevector.h"
#include "paths/KmerPath.h"

/// Generate paths given a working directory and a list fastb files.
/// readsInHead is the list of fastb filenames but without the .fastb extension.
/// The resulting paths will have the a .paths.kN extension.
void CommonPather(const int K, 
                  const String& dirIn,
		  const vec<String>& readsInHead,
                  const int NUM_THREADS,
                  const String CHECKPOINT_HEAD = "");

/// Generate paths given a working directory and a list fastb files and paths filenames.
/// readsIn is the list of fastb files to path.
/// pathsOut is the list of paths filenames (to which .kN will be appended).
void CommonPather(const int K, 
                  const String& dirIn, 
                  const vec<String>& readsIn,
		  const vec<String>& pathsOut,
                  const int NUM_THREADS,
                  const String CHECKPOINT_HEAD = "");

/// Generate paths given a list fastb files and path filenames.
/// readsIn is the list of fastb files to path.
/// pathsOut is the list of paths filenames (to which .kN will be appended).
void CommonPather(const int K, 
                  const vec<String>& readsIn,
                  const vec<String>& pathsOut, 
                  const int NUM_THREADS,
                  const String CHECKPOINT_HEAD = "",
                  const Bool appendKSize = True, 
                  const Bool verify = False);

/// Generate paths from a set of vecbasevectors - no external files required
/// readsIn is the vec of pointers to the vecbasevectors to path.
/// pathsOut is the vec of pointers to the veckmerpaths in which to store the paths
void CommonPather(const int K, 
                  const vec<const vecbvec*>& readsIn, 
                  const vec<vecKmerPath*>& pathsOut,
                  const int NUM_THREADS,
                  const String CHECKPOINT_HEAD = "");

/// Wrapper that deals with vec of const/non-const pointer issues.
void CommonPather(const int K, 
                  const vec<vecbvec*>& readsIn, 
                  const vec<vecKmerPath*>& pathsOut,
                  const int NUM_THREADS,
                  const String CHECKPOINT_HEAD = "");

/// Convenience functions to path two or three vecbasevectors together:

/// Generate paths given two vecbasevectors
/// readsIn1, readsIn2 are the vecbasevectors to path.
/// pathsOut1, pathsOut2 are the veckmerpaths in which to store the paths
void CommonPather(const int K, 
                  const vecbvec& readsIn1, 
                  const vecbvec& readsIn2,
                  vecKmerPath& pathsOut1, 
                  vecKmerPath& pathsOut2,
                  const int NUM_THREADS,
                  const String CHECKPOINT_HEAD = "");

/// Generate paths given three vecbasevectors
/// readsIn1, readsIn2, readsIn3 are the vecbasevectors to path.
/// pathsOut1, pathsOut2, pathsOut3 are the veckmerpaths in which to store the paths
void CommonPather(const int K,
                  const vecbvec& readsIn1, 
                  const vecbvec& readsIn2, 
                  const vecbvec& readsIn3,
                  vecKmerPath& pathsOut1, 
                  vecKmerPath& pathsOut2, 
                  vecKmerPath& pathsOut3,
                  const int NUM_THREADS,
                  const String CHECKPOINT_HEAD = "");


#endif
