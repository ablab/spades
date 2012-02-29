/**
   \file

   Extraction and sorting of kmers; see SortKmers().  This file exists just for explicit template instantiation; the actual
   kmer extraction and sorting code is in SortKmersImpl.h .

   \ingroup grp_kmerGathering
*/


#include "kmers/SortKmersImpl.h"
#include "kmers/KmerShape.h"

// Include files to link in instantiations.

#include "sort_kmers/SortKmersA.h"
#include "sort_kmers/SortKmersB.h"
#include "sort_kmers/SortKmersC.h"
#include "sort_kmers/SortKmersD.h"

void DEFAULT_SORT_KMERS_END_PASS( int pass ) {
  Dot( cout, pass );
}

void (*SORT_KMERS_END_PASS)( int pass ) = DEFAULT_SORT_KMERS_END_PASS;




