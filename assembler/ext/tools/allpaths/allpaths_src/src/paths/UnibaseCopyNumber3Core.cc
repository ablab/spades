///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Note: this code contains the guts of the copy number computations.  There is
// also code in UnipathPatcher.cc that follows essentially the same method.
// Thus if the code for the main routine UnibaseCopyNumber3Core here is changed,
// we should check to see if UnipathPatcher.cc needs to be changed as well.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "Basevector.h"
#include "CoreTools.h"
#include "kmers/KmerParcelsTools.h"
#include "paths/PdfEntry.h"
#include "paths/UnibaseCopyNumberCommon.h"
#include "paths/UnibaseCopyNumber3Core.h"

void ComputeProbs(
     // INPUTS:
     const int K,
     const int PLOIDY,
     double occCnPloidy,
     double THRESH,
     double ERR_RATE,
     const vecbasevector& unibases,
     const vec<int64_t>& n_kmer_hits,
     // OUTPUTS:
     VecPdfEntryVec& cn_pdfs,
     vec<double>& cn_raw,
     vec<longlong>& cnHisto,
     // LOGGING:
     const Bool VERBOSE,
     ostream& logout )
{
  // Normalize the total unibase coverage into an average.

  logout << Date() << ": Computing averages" << endl;
  vec<double> avg_kmer_freq;
  for ( size_t i = 0; i < unibases.size(); i++ ) {
    double n_kmers = unibases[i].size() - K + 1;
    avg_kmer_freq.push_back( n_kmer_hits[i] / n_kmers );
  }

  cn_pdfs.clear( );
  cn_pdfs.resize( unibases.size( ) );

  cout << Date( ) << ": starting loop through unibases" << endl;
  #pragma omp parallel for
  for ( size_t u = 0; u < unibases.size(); ++u ) {
    PdfEntryVec cprob, copyno_prob;
    longlong n =  (int) round( (double)avg_kmer_freq[u] );
    int nKmers = unibases[u].isize() -K +1;
    ForceAssertGe( n, 0 );
    vec<double> log_nfact( 1, 0 );
    int lnfacts = log_nfact.size( );
    if ( n + 1 > lnfacts ){    
      log_nfact.resize( n + 1 );
      for ( int i = lnfacts; i <= n; i++ )
	log_nfact[i] = log_nfact[i-1] + log( double(i) );    
    }
    double q          = occCnPloidy / (double) PLOIDY;
    double thresh     = THRESH;
    double error_rate = ERR_RATE;
    double logplast   = 0.0, sump = 0.0;
    double logthresh  = log(thresh);
    double logpmax    = -numeric_limits<double>::max();
    
    cn_raw[u] = (double)avg_kmer_freq[u] / q;
    /// This part is a holdover from a previous code. Technically it is not correct
    /// since we are using an average of a sample now.

    double error_prior = 1 - pow( 1-error_rate, K );  
    // naive probability that the a K-mer contains an error;
    // used as a prior against the C=0 hypothesis.
    
    // First model C=0 (ie error kmer):
    if( ERR_RATE > 0 ) {
      double c = ERR_RATE / 3;
      //double qm = q + 5 * sqrt(q);
      double logp = -q * c + double( n >= 1 ? n-1 : 0 ) * log( q * c ) - log_nfact[n >= 1 ? n-1 : 0];
      // use n-1 instead of n, since clearly it does appear once
      double p = error_prior * exp(logp);
      cprob.push_back( pdf_entry( 0, logp ) );    
    }
    // Now model any positive copy number; prior is all are equally likely.
    for ( int C = 1; ; C++ ){    
      if ( (C % PLOIDY) != 0 && C > 2 ) continue;
      double c = C;
      double logp = -q * c + double(n) * log( q * c ) - log_nfact[n];
      if ( (C>PLOIDY && logp < logplast) && logp <= logpmax + logthresh ) break;
      cprob.push_back( pdf_entry( C, logp ) );
      logplast = logp;
      logpmax = Max( logpmax, logp );
    }
    // Let prob(C=0) into calculation of pmax, if needed; harmless if not
    logpmax = Max( logpmax, cprob[0].second );
    
    for ( PdfEntryVec::size_type i = 0; i < cprob.size( ); i++ )
      if ( cprob[i].second >= logpmax + logthresh ) copyno_prob.push_back( cprob[i] );
    
    double pmax1 = cprob.at(0).second;
    double cnmax1 = cprob[0].first;
    for ( PdfEntryVec::size_type i = 1; i < cprob.size( ); i++ )
      if ( cprob[i].second > pmax1 ){
	pmax1 = cprob[i].second;
	cnmax1 = cprob[i].first;
      }
 
    // Approximate normalization in log space, using logpmax.

    for ( PdfEntryVec::size_type i = 0; i < copyno_prob.size( ); i++ )
      copyno_prob[i].second = exp(copyno_prob[i].second - logpmax);
    for ( PdfEntryVec::size_type i = 0; i < copyno_prob.size( ); i++ )
      sump += copyno_prob[i].second;
    for ( PdfEntryVec::size_type i = 0; i < copyno_prob.size( ); i++ )
      copyno_prob[i].second /= sump;    

    if (VERBOSE) {
      double maxProb = -1.0;
      int maxProbCN = -1;
      for ( unsigned i = 0; i < copyno_prob.size(); ++i )
	if ( copyno_prob[i].Prob() > maxProb ) {
	  maxProb   = copyno_prob[i].Prob();
	  maxProbCN = copyno_prob[i].NumCopies();
	}
      #pragma omp critical
      { 
	if ( maxProbCN > cnHisto.isize() -1 )
	  cnHisto.resize( maxProbCN + 1 );
      }
      cnHisto[ maxProbCN]++;
    }
    
    cn_pdfs[u] = copyno_prob;
  }

}

void UnibaseCopyNumber3Core( 
     // INPUTS:
     const int K, const vecbasevector& reads, const vecbasevector& unibases, 
     const int NUM_THREADS, const String& reads_head, const String& unibases_head,
     const int PLOIDY,
     // LOGGING CONTROL:
     const Bool VERBOSE, ostream& logout,
     // OUTPUT:
     VecPdfEntryVec& cn_pdfs,
     vec<double>& cn_raw,
     vec<int64_t>& n_kmer_hits_raw,
     // HEURISTIC PARAMETERS:
     const double THRESH, const double ERR_RATE, const Bool CORRECT_BIAS )
{
  // Make coupled KmerParcels out of the reads and unibases.

  logout << Date() << ": Creating coupled KmerParcels out of the reads and unibases" << endl;


  // ---- build coupled parcels for reads and unibases
  KmerParcelsDiskStore parcels_reads(K, reads_head);
  KmerParcelsDiskStore parcels_unibases(K, unibases_head);

  KmerParcelsCoupledBuilder builder(reads, unibases, parcels_reads, parcels_unibases,
                                    NUM_THREADS);
  builder.Build();
  size_t num_parcels = parcels_reads.GetNumParcels();
  

  /// get a small number of longest unibases and estimate the expected occurrence of kmers

  double occCnPloidy = -1.0;
  vec<longlong> biasOccs( K+1, 1), biasInst( K+1, 0);
  {

    vec<unsigned> ulens( unibases.size(), 0 );
    for ( size_t i = 0; i < unibases.size( ); i++ )
      ulens[ i ] = unibases[ i ].size();

    int min_length, nlongest;
    GetMinLength( unibases, min_length, nlongest );

    // Open up the coupled KmerParcels and process the batches,
    // finding the (shared) basevector and the frequency in each dataset.
    // We want to find the total number of times each kmer from each unibase
    // appears in the reads.

    logout << Date() << ": Preliminary bias computation - Parsing the coupled KmerParcels" << endl;
    n_kmer_hits_raw.resize_and_set( unibases.size(), 0 );
    vec<int64_t> n_kmer_hits( unibases.size( ), 0 );
    
    BaseVec kmer_bv;
    for (size_t parcel_ID = 0; parcel_ID != num_parcels; parcel_ID++) {
      KmerParcelReader parcel_reads   (parcels_reads,    parcel_ID);
      KmerParcelReader parcel_unibases(parcels_unibases, parcel_ID);
      
      while (GetNextMatchingKmerBatches(parcel_reads, parcel_unibases)) {
        const KmerBatch & batch_reads    = parcel_reads.CurrentKmerBatch();
        const KmerBatch & batch_unibases = parcel_unibases.CurrentKmerBatch();
	
	// For each time this kmer appears on a unibase, add the kmer's
	// read frequency to that unibase's coverage.

	size_t kmer_freq_reads = batch_reads.GetKmerFreq();
        batch_reads.GetKmerBaseVec(&kmer_bv);
	int gc = kmer_bv.GcBases();
	
	for ( size_t i = 0; i < batch_unibases.GetKmerFreq(); i++ ) {
	  size_t uID = batch_unibases[i].GetReadID();
          n_kmer_hits_raw[uID] += kmer_freq_reads;
	  // Take only the 100 longest unipaths.
	  if ( (int) ulens[uID] >= min_length ) {
            n_kmer_hits[uID] += kmer_freq_reads;
            biasOccs[gc] += kmer_freq_reads;
            biasInst[gc]++;
          }
	}
      }
    }
    
    vec<double> avg_kmer_freq;
    for ( size_t i = 0; i < unibases.size(); i++ ) {
      double n_kmers = unibases[i].size() - K + 1;
      avg_kmer_freq.push_back( n_kmer_hits[i] / n_kmers );
    }
    
    occCnPloidy = double(Sum(avg_kmer_freq)) / (double) nlongest; 
  }








  
  vec<double> biasCurveLoc( K+1, 1.0);
  if ( CORRECT_BIAS )
    ComputeBias( K, occCnPloidy, biasOccs, biasInst, biasCurveLoc, VERBOSE, logout );

  // Open up the coupled KmerParcels and process the batches,
  // finding the (shared) basevector and the frequency in each dataset.
  // We want to find the total number of times each kmer from each unibase
  // appears in the reads - this time, adjusted for bias!

  logout << Date() << ": Parsing the coupled KmerParcels" << endl;
  vec<int64_t> n_kmer_hits( unibases.size(), 0 );
  
  BaseVec kmer_bv;
  for (size_t parcel_ID = 0; parcel_ID != num_parcels; parcel_ID++ ) {
    KmerParcelReader parcel_reads   (parcels_reads,    parcel_ID);
    KmerParcelReader parcel_unibases(parcels_unibases, parcel_ID);
    
    while (GetNextMatchingKmerBatches(parcel_reads, parcel_unibases)) {
      const KmerBatch & batch_reads    = parcel_reads.CurrentKmerBatch();
      const KmerBatch & batch_unibases = parcel_unibases.CurrentKmerBatch();
      
      // For each time this kmer appears on a unibase, add the kmer's
      // read frequency (adjusted for bias) to that unibase's coverage.
      
      batch_reads.GetKmerBaseVec(&kmer_bv);
      int gc = kmer_bv.GcBases();
	
      double freq = (double)batch_reads.GetKmerFreq() / biasCurveLoc[ gc ];
      
      for ( size_t i = 0; i < batch_unibases.GetKmerFreq(); i++ )
	n_kmer_hits[ batch_unibases[i].GetReadID() ] += freq;
    }
  }

  vec<longlong> cnHisto;

  ComputeProbs( K, PLOIDY, occCnPloidy, THRESH, ERR_RATE, unibases, n_kmer_hits,
		cn_pdfs, cn_raw, cnHisto, VERBOSE, logout );

  if ( VERBOSE ) {
    logout << "copy number histogram\n";
    cnHisto.Println(logout);
  }

}
