/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "feudal/BaseVec.h"
#include "kmers/KmerRecord.h"
#include "VecTemplate.h"

int const kmer_with_count_base::max_count;

FOR_ALL_K(INSTANTIATE_KMER_RECORD_FOR_K,unused);

template<int K> void kmer<K>::SetToSubOf( const basevector& source, 
     const size_type start )
{   size_t len = K;
    AssertLe( start, source.size() );
    AssertLe( len, source.size()-start );
    size_t end = len & ~15;
    int32_t* dst = (int32_t*) &data_;
    for ( size_t idx = 0; idx < end; idx += 16)
         *dst++ = source.extractKmer(start+idx, 16);
    if ( end < len ) *dst = source.extractKmer(start+end, len-end);   }

template void kmer<4>::SetToSubOf( const basevector&, const size_type );
template void kmer<20>::SetToSubOf( const basevector&, const size_type );
template void kmer<40>::SetToSubOf( const basevector&, const size_type );
template void kmer<80>::SetToSubOf( const basevector&, const size_type );
template void kmer<96>::SetToSubOf( const basevector&, const size_type );
template void kmer<100>::SetToSubOf( const basevector&, const size_type );
template void kmer<400>::SetToSubOf( const basevector&, const size_type );
template void kmer<640>::SetToSubOf( const basevector&, const size_type );
