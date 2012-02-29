///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef ASSEMBLY_CLEANUP_TOOLS_H
#define ASSEMBLY_CLEANUP_TOOLS_H

#include "CoreTools.h"
#include "efasta/EfastaTools.h"
#include "VecUtilities.h"

/**
 * Class: Assembly
 *
 * This is a class for manimpulating a superb and associate contigs file).
 */

class Assembly{
  
 public:
  
  vec<superb> scaffolds;
  vec<efasta> efastas;

  vec<String> scaffMap;
  vec<String> tigMap;
  
  // -------------- Constructors
  Assembly( const vec<superb>& scaffolds, const vec<efasta>& efastas );
  Assembly( const String scaffoldsFile, const String efastaFile );

  // --- useful functions
  void check_integrity() const;
  size_t scaffoldsTotLen() const;
  size_t scaffoldsRedLen() const;
  size_t scaffoldsNtigs() const;

  void remove_small_scaffolds( const int MIN_SCAFFOLD_SIZE );
  void remove_contigs( const vec<Bool>& to_remove );
  void remove_small_contigs( const int MIN_CONTIG_SIZE_SOLO, const int MIN_CONTIG_SIZE_IN );
  void remove_unused_contigs();
  void dedup();
  void dedup2();
  void reorder();
  // renumber all the contigs sequentially according to the scaffold
  void renumber();

  vec<String> getScaffMap() const { return scaffMap; }
  vec<String> getTigMap() const { return tigMap; }

  void Write( const String head_out ) const;
  void WriteExtra( const String head_out ) const;
  void WriteAll( const String head_out ) const {
    Write( head_out );
    WriteExtra( head_out );
  }

};

#endif
