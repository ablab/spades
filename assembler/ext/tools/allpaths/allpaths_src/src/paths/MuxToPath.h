// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

#ifndef PATHS_MUXTOPATH_H
#define PATHS_MUXTOPATH_H

#include "paths/Mux.h"
#include "paths/KmerPath.h"

/// A class which encapsulates the conversion of a vec<Mux> to
/// its corresponding kmer path.

class MuxToPath {
public:
  MuxToPath( const vecKmerPath* pathsFw, const vecKmerPath* pathsRc ) :
    mp_pathsFw( pathsFw ), mp_pathsRc( pathsRc ) { }
  // It should be OK to instantiate this class with garbage pointers, 
  // as long as you never call any member functions.

  // If the kmer_path is empty, then the read corresponding to the 
  // first Mux in the vec is copied over in its entirety.
  // Otherwise, the first mux should extend the end of kmer_path.
  void ToPath( const vec<Mux>& muxes, 
	       const KmerPath& path_so_far, 
	       KmerPath& ans ) const;

  void ExtendByMux( const Mux& the_mux, 
		    const KmerPath& path_so_far,
		    KmerPath& ans ) const;

  bool ExtendByKmersIfMatch( const OrientedKmerPathId& okpid,
			     const int numKmers,
			     const KmerPath& path_so_far,
			     KmerPath& ans ) const;

private:
  const vecKmerPath* mp_pathsFw;
  const vecKmerPath* mp_pathsRc;

};

#endif
