///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////



#include "Basevector.h"
#include "Qualvector.h"
#include "reporting/PerfStat.h"



static inline 
String Tag(String S = "KSC") { return Date() + " (" + S + "): "; } 



#include "kmers/KmerSpectra.h"
#include "kmers/naif_kmer/NaifKmerizer.h" 
#include "kmers/naif_kmer/KernelKmerSpectralizer.h" 






void kmer_spectrum_compute(const BaseVecVec & bvv,
                           KmerSpectrum     * p_kspec,
                           const unsigned     VERBOSITY, 
                           const unsigned     NUM_THREADS,
                           const size_t       mem_mean_ceil)
{
  const unsigned K = p_kspec->K();

  if (K <= 29) {
    KernelKmerSpectralizer<Kmer29> spectralizer(bvv, K, p_kspec);
    naif_kmerize(&spectralizer, NUM_THREADS, VERBOSITY, mem_mean_ceil);
  }
  else if (K <= 60) {
    KernelKmerSpectralizer<Kmer60> spectralizer(bvv, K, p_kspec);
    naif_kmerize(&spectralizer, NUM_THREADS, VERBOSITY, mem_mean_ceil);
  }
  else if (K <= 124) {
    KernelKmerSpectralizer<Kmer124> spectralizer(bvv, K, p_kspec);
    naif_kmerize(&spectralizer, NUM_THREADS, VERBOSITY, mem_mean_ceil);
  }
  else if (K <= 248) {
    KernelKmerSpectralizer<Kmer248> spectralizer(bvv, K, p_kspec);
    naif_kmerize(&spectralizer, NUM_THREADS, VERBOSITY, mem_mean_ceil);
  }
  else if (K <= 504) {
    KernelKmerSpectralizer<Kmer504> spectralizer(bvv, K, p_kspec);
    naif_kmerize(&spectralizer, NUM_THREADS, VERBOSITY, mem_mean_ceil);
  }
  else {
    cout << Tag() << "K= " << K << " not implemented." << endl;
    exit(1);
  }
   
        
}
    








void kmer_spectra_by_affixes_compute(const BaseVecVec   & bvv,
                                     KmerSpectraAffixes * p_kspecs,
                                     const Validator    & validator,
                                     const unsigned       VERBOSITY,
                                     const unsigned       NUM_THREADS,
                                     const size_t         mem_mean_ceil)
{
  const unsigned K = p_kspecs->K();

  if (K <= 29) {
    KernelKmerAffixesSpectralizer<Kmer29> spectralizer(bvv, K, p_kspecs, &validator);
    naif_kmerize(&spectralizer, NUM_THREADS, VERBOSITY, mem_mean_ceil);
  }
  else if (K <= 60) {
    KernelKmerAffixesSpectralizer<Kmer60> spectralizer(bvv, K, p_kspecs, &validator);
    naif_kmerize(&spectralizer, NUM_THREADS, VERBOSITY, mem_mean_ceil);
  }
  else if (K <= 124) {
    KernelKmerAffixesSpectralizer<Kmer124> spectralizer(bvv, K, p_kspecs, &validator);
    naif_kmerize(&spectralizer, NUM_THREADS, VERBOSITY, mem_mean_ceil);
  }
  else if (K <= 248) {
    KernelKmerAffixesSpectralizer<Kmer248> spectralizer(bvv, K, p_kspecs, &validator);
    naif_kmerize(&spectralizer, NUM_THREADS, VERBOSITY, mem_mean_ceil);
  }
  else if (K <= 504) {
    KernelKmerAffixesSpectralizer<Kmer504> spectralizer(bvv, K, p_kspecs, &validator);
    naif_kmerize(&spectralizer, NUM_THREADS, VERBOSITY, mem_mean_ceil);
  }
  else {
    cout << Tag() << "K= " << K << " not implemented." << endl;
    exit(1);
  }
}





void kmer_bi_spectrum_compute(const BaseVecVec   & bvv,
                              const size_t         nbv1,
                              KmerBiSpectrum     * p_kbspec,
                              const unsigned       VERBOSITY,
                              const unsigned       NUM_THREADS,
                              const size_t         mem_mean_ceil)
{
  const unsigned K = p_kbspec->K();

  if (K <= 29) {
    KernelKmerBiSpectralizer<Kmer29> spectralizer(bvv, nbv1, K, p_kbspec);
    naif_kmerize(&spectralizer, NUM_THREADS, VERBOSITY, mem_mean_ceil);
  }
  else if (K <= 60) {
    KernelKmerBiSpectralizer<Kmer60> spectralizer(bvv, nbv1, K, p_kbspec);
    naif_kmerize(&spectralizer, NUM_THREADS, VERBOSITY, mem_mean_ceil);
  }
  else if (K <= 124) {
    KernelKmerBiSpectralizer<Kmer124> spectralizer(bvv, nbv1, K, p_kbspec);
    naif_kmerize(&spectralizer, NUM_THREADS, VERBOSITY, mem_mean_ceil);
  }
  else if (K <= 248) {
    KernelKmerBiSpectralizer<Kmer248> spectralizer(bvv, nbv1, K, p_kbspec);
    naif_kmerize(&spectralizer, NUM_THREADS, VERBOSITY, mem_mean_ceil);
  }
  else if (K <= 504) {
    KernelKmerBiSpectralizer<Kmer504> spectralizer(bvv, nbv1, K, p_kbspec);
    naif_kmerize(&spectralizer, NUM_THREADS, VERBOSITY, mem_mean_ceil);
  }
  else {
    cout << Tag() << "K= " << K << " not implemented." << endl;
    exit(1);
  }
        
}






void kmer_spectra_by_quals_compute(const BaseVecVec   & bvv,
				   const QualVecVec   & qvv,
				   KmerQualitySpectra * p_kqspec,
				   const unsigned       VERBOSITY,
				   const unsigned       NUM_THREADS,
                                   const size_t         mem_mean_ceil)
{
  const unsigned K = p_kqspec->K();

  if (K <= 29) {
    KernelKmerQualitySpectralizer<Kmer29> spectralizer(bvv, qvv, K, p_kqspec);
    naif_kmerize(&spectralizer, NUM_THREADS, VERBOSITY, mem_mean_ceil);
  }
  else if (K <= 60) {
    KernelKmerQualitySpectralizer<Kmer60> spectralizer(bvv, qvv, K, p_kqspec);
    naif_kmerize(&spectralizer, NUM_THREADS, VERBOSITY, mem_mean_ceil);
  }
  else if (K <= 124) {
    KernelKmerQualitySpectralizer<Kmer124> spectralizer(bvv, qvv, K, p_kqspec);
    naif_kmerize(&spectralizer, NUM_THREADS, VERBOSITY, mem_mean_ceil);
  }
  else if (K <= 248) {
    KernelKmerQualitySpectralizer<Kmer248> spectralizer(bvv, qvv, K, p_kqspec);
    naif_kmerize(&spectralizer, NUM_THREADS, VERBOSITY, mem_mean_ceil);
  }
  else if (K <= 504) {
    KernelKmerQualitySpectralizer<Kmer504> spectralizer(bvv, qvv, K, p_kqspec);
    naif_kmerize(&spectralizer, NUM_THREADS, VERBOSITY, mem_mean_ceil);
  }
  else {
    cout << Tag() << "K= " << K << " not implemented." << endl;
    exit(1);
  }
}








void genome_analysis_report(const KmerSpectrum & kspec,
                            const unsigned       read_len,
                            const unsigned       PLOIDY,
                            const unsigned       KF_LOW,
                            const unsigned       VERBOSITY)
{
  cout << Tag() << "Estimating genome size." << endl; 
  kspec.analyze(PLOIDY, read_len, KF_LOW, VERBOSITY);

  const size_t G = kspec.genome_size();
  const bool failed = (G == 0);

  const size_t G1      = (failed) ? 0 : kspec.genome_size_unique();
  const size_t GR      = (failed) ? 0 : kspec.genome_size_repetitive();
  const float coverage = (failed) ? 0 : kspec.coverage();

  cout << Tag() << "------------------- Kmer Spectrum Analysis -------------------" << endl;
  cout << Tag() << "Genome size estimate        = " 
       << setw(14) << ToStringAddCommas(G) << " bases" << endl;
  cout << Tag() << "Genome size estimate CN = 1 = " 
       << setw(14) << ToStringAddCommas(G1) << " bases ( " 
       << fixed << setw(5) << setprecision(1) 
       << (failed ? 0 : 100 * float(G1)/float(G)) << " % )" << endl;
  cout << Tag() << "Genome size estimate CN > 1 = " 
       << setw(14) << ToStringAddCommas(GR) << " bases ( "
       << fixed << setw(5) << setprecision(1) 
       << (failed ? 0 : 100 * float(GR)/float(G)) << " % )" << endl;
        
  // ---- computing coverage
        
  cout << Tag() << "Coverage estimate           = " 
       << setw(14) << setprecision(0) << coverage << " x" << endl;


  // ---- bias standard deviation (coverage independent)

  const float sd_bias = failed ? 0 : kspec.bias_stddev();
  cout << Tag() << "Bias stddev at scale > K    = " 
       << setw(14) << setprecision(2) << sd_bias << endl;


  // ---- estimating base rate based on number of bad kmers

  const float pe = kspec.fraction_of_error_kmers() / kspec.K();
  const float Q = -10.0 * log(pe) / log(10);
  cout << Tag() << "Base error rate estimate    = " 
       << setw(14) << setprecision(4) << pe 
       << " (Q = " << setprecision(1) << Q << ")" << endl;
    

        
  // ---- ploidy report

  if (PLOIDY == 2) {
    cout << Tag() << "Ploidy                      = " << setw(14) << PLOIDY << endl; 
    const size_t d_SNPs = failed ? 0 : kspec.d_SNP();
    cout << Tag() << "SNP rate                   >= " 
	 << setw(14) << ("1/" + ToString(d_SNPs)) << endl;

    const float p = failed ? 0.0 : 1.0 / d_SNPs;
    const float p_close = 1.0 - pow(1.0 - p, static_cast<int>(kspec.K()));
    cout << Tag() << "SNPs closer than K         >= " 
         << setw(14) << setprecision(0) << 100 * p_close * (2 - p_close) << " %" << endl;
  }
  else {
    cout << Tag() << "SNP rate not computed (PLOIDY = " << PLOIDY << ")." << endl;
  } 
  cout << Tag() << "--------------------------------------------------------------" << endl;

}







void genome_size_perf_stats(const KmerSpectrum & kspec,
                            const unsigned       read_len,
                            const unsigned       PLOIDY,
                            const unsigned       KF_LOW,
                            const unsigned       VERBOSITY)
{
  kspec.analyze(PLOIDY, read_len, KF_LOW, VERBOSITY);
  
  const size_t G = kspec.genome_size();
  const bool failed = (G == 0);

  const size_t GR      = (failed) ? 0 : kspec.genome_size_repetitive();
  const float coverage = (failed) ? 0 : kspec.coverage();
  const float sd_bias  = failed ? 0 : kspec.bias_stddev();

  PerfStat::log() << std::fixed << std::setprecision(0)
                  << PerfStat( "genome_size_est", 
                               "estimated genome size in bases", 
                               G);
  PerfStat::log() << std::fixed << std::setprecision(1)
                  << PerfStat( "genome_repetitiveness_est", 
                               "% genome estimated to be repetitive (at K=" + ToString(kspec.K()) + " scale)", 
                               (failed ? 0 : 100 * GR / G));
  PerfStat::log() << std::fixed << std::setprecision(0)
                  << PerfStat( "genome_cov_est", 
                               "estimated genome coverage by fragment reads", 
                               coverage);
  PerfStat::log() << std::fixed << std::setprecision(2)
                  << PerfStat( "bias_stddev_est", 
                               "estimated standard deviation of sequencing bias (at K=" + ToString(kspec.K()) + " scale)", 
                               sd_bias);

}










