///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
//  Author: Filipe Ribeiro - Nov 2011 - <crdhelp@broadinstitute.org>
//

#ifndef KMERS__KMER_SPECTRUM_CORE_H
#define KMERS__KMER_SPECTRUM_CORE_H


#include "Basevector.h"
#include "Qualvector.h"
#include "kmers/KmerSpectra.h"

#include "kmers/naif_kmer/Kmers.h"



void kmer_spectrum_compute(const BaseVecVec & bvv,
                           KmerSpectrum     * p_kspec,
                           const unsigned     VERBOSITY, 
                           const unsigned     NUM_THREADS,
                           const size_t       mem_mean_ceil = 0);

void kmer_spectra_by_affixes_compute(const BaseVecVec   & bvv,
                                     KmerSpectraAffixes * p_kspecs,
                                     const Validator    & validator,
                                     const unsigned       VERBOSITY,
                                     const unsigned       NUM_THREADS,
                                     const size_t         mem_mean_ceil = 0);

void kmer_bi_spectrum_compute(const BaseVecVec   & bvv,
                              const size_t         nbv1,
                              KmerBiSpectrum     * p_kbspec,
                              const unsigned       VERBOSITY,
                              const unsigned       NUM_THREADS,
                              const size_t         mem_mean_ceil = 0);


void kmer_spectra_by_quals_compute(const BaseVecVec   & bvv,
				   const QualVecVec   & qvv,
				   KmerQualitySpectra * p_kqspec,
				   const unsigned       VERBOSITY,
				   const unsigned       NUM_THREADS,
                                   const size_t         mem_mean_ceil = 0);

void genome_analysis_report(const KmerSpectrum & kspec,
                            const unsigned       read_len,
                            const unsigned       PLOIDY,
                            const unsigned       KF_LOW,
                            const unsigned       VERBOSITY);

void genome_size_perf_stats(const KmerSpectrum & kspec,
                            const unsigned       read_len,
                            const unsigned       PLOIDY,
                            const unsigned       KF_LOW,
                            const unsigned       VERBOSITY);

#endif
