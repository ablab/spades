// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

#include "paths/KmerPathMuxSearcher.h"


// This just calls the templatized SearchDirector with the
// appropriate SearchAgent
void KmerPathMuxSearcher::FindClosures( const int closingReadId,
					const int openingReadId, 
					const int minAcceptableExtLength,
					const int maxAcceptableExtLength,
					MuxSearchResult &result ) const
{
  if( m_policies.empty() && m_min_perfect_match <= 1 ) {
    MuxSearchAgentSimple msa( mp_muxGraph, mp_readfillDB, mp_subList, m_verbosity );
    SearchDirector( msa, openingReadId, closingReadId,
		    minAcceptableExtLength, maxAcceptableExtLength,
		    result );
  }
  else {
    MuxSearchAgentGoodReads msa( mp_muxGraph, mp_readfillDB, 
				 mp_subList, mp_readLengths, m_max_read_length,
				 m_min_perfect_match,
				 m_verbosity, m_policies );
    SearchDirector( msa, openingReadId, closingReadId,
		    minAcceptableExtLength, maxAcceptableExtLength,
		    result );
  }
}




