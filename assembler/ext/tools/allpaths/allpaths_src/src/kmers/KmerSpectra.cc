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

#include "Basevector.h"
#include "kmers/KmerSpectra.h"


size_t KmerSpectrum::sum() const
{
  const vec<size_t> & ndk = *this;     // number of distinct kmers at freq kf
  const size_t nkf = ndk.size();
  size_t sum = 0;
  for (size_t kf = 0; kf != nkf; kf++)
    sum += ndk[kf];
  return sum;  
}





void KmerSpectrum::analyze(const unsigned ploidy, 
                           const unsigned read_len,
                           const size_t kf_min1_arg,
                           const unsigned verbosity) const 
{
  // ---- reset genome statistics

  _genome_size = 0;
  _genome_size_unique = 0;
  _genome_size_repetitive = 0;
  _d_SNP = 0;
  _stddev_bias = 0;


  const vec<size_t> & ndk = *this;     // number of distinct kmers at freq kf
  const size_t nkf = ndk.size();

  vec<size_t> cndk(nkf, 0);   // cumulative number of distinct kmers at freq kf
    
  vec<size_t> nk(nkf, 0);    // number of kmers at freq kf
  vec<size_t> cnk(nkf, 0);    // cumulative number of kmers at freq kf

  const size_t kf_ceil = nkf - 1;
  if (verbosity > 1) cout << "kf_ceil         = " << kf_ceil << endl;


  for (size_t kf = 0; kf != nkf; kf++)
    nk[kf] = kf * ndk[kf];
    
  for (size_t kf = 1; kf != nkf; kf++) {
    cndk[kf] = cndk[kf - 1] + 0.5 * (ndk[kf - 1] + ndk[kf]);
    cnk [kf] = cnk [kf - 1] + 0.5 * (nk [kf - 1] + nk [kf]);
  }


  // ---- Separate kmer spectrum in 5 regions based on the kf:
  //    1       ... kf_min1   : bad kmers with low frequency
  //    kf_min1 ... kf_min2   : good kmers CN = 1/2 (SNPs)
  //    kf_min2 ... kf_min3   : good kmers CN = 1 
  //    kf_min3 ... kf_hi     : good kmers CN > 1 (repetitive)
  //    kf_hi   ... inf       : bad kmers with high frequency
 
 
  // ---- min1: find first minimum
  
  if (kf_min1_arg > kf_ceil) {
    cout << "**** Kmer frequencies are too low. Bad data?" << endl;
    return;
  }

  _kf_min1 = kf_min1_arg;
  // first search lower frequencies
  while (_kf_min1 - 1 >= 2 && nk[_kf_min1 - 1] < nk[_kf_min1]) 
    _kf_min1--;
  // next search higher frequencies
  while (_kf_min1 <= kf_ceil && nk[_kf_min1 + 1] < nk[_kf_min1]) 
    _kf_min1++;

  if (_kf_min1 > kf_ceil) {
    cout << "**** Could not find kmer spectrum local minimum. Bad data?" << endl;
    return;
  }
  
  if (verbosity > 1) cout << "kf_min1         = " << _kf_min1 << endl;


  // ---- max2: find absolute maximum max2 above first minimum min1
    
  _kf_max2 = _kf_min1;
  for (size_t kf = _kf_min1 + 1; kf < 0.8 * kf_ceil ; kf++)
    if (nk[kf] > nk[_kf_max2])
      _kf_max2 = kf;

  if (verbosity > 1) cout << "kf_max2         = " << _kf_max2 << endl;

  if (_kf_max2 < kf_min1_arg) {
    cout << "**** Coverage too low. Can't estimate genome size." << endl;
    return;
  }
    
  // ---- max2: resetting max2 for cases of very high polymorphism
  if (ploidy == 2) {
    size_t ndk_half   = ndk[_kf_max2 / 2];
    size_t ndk_double = ndk[_kf_max2 * 2];
    if (ndk_double > ndk_half) _kf_max2 *= 2;
      
    if (verbosity > 1) cout << "kf_max2         = " << _kf_max2 << endl;
  }


  // ---- max1: SNPs local maximum max1 as half global maximum max2

  _kf_max1 = _kf_max2 / 2;

  if (verbosity > 1) cout << "kf_max1         = " << _kf_max1 << endl;
    
  // ---- min2: SNPs local minimum min2 between max1 and max2

  _kf_min2 = _kf_max1 * (2 * ndk[_kf_max1] + ndk[_kf_max2]) / (ndk[_kf_max1] + ndk[_kf_max2]);

  if (verbosity > 1) cout << "kf_min2         = " << _kf_min2 << endl;


  // ---- min1: refine between min1 and max2/2

  for (size_t kf = _kf_min1 + 1; kf < _kf_max1; kf++)
    if (nk[kf] < nk[_kf_min1])
      _kf_min1 = kf;


    
  // ---- haploid? 
  if (ploidy != 2) {
    _kf_max1 = _kf_min2 = _kf_min1;
      
    if (verbosity > 1) cout << "kf_max1         = " << _kf_max1 << endl;
    if (verbosity > 1) cout << "kf_min2         = " << _kf_min2 << endl;
  }



  // ---- min3: not a minimum, really. upper edge of main peak
    
  _kf_min3 = _kf_max2 * 3 / 2;
  if (_kf_min3 > kf_ceil) {
    cout << "**** Error analysing spectrum: kf_min3 (" << _kf_min3 << ") > kf_ceil (" << kf_ceil << ")." << endl;
    return;
  }

  if (verbosity > 1) cout << "kf_min3         = " << _kf_min3 << endl;

  if (verbosity > 1) cout << "kfs             :" 
                    << " " << _kf_min1
                    << " " << _kf_max1
                    << " " << _kf_min2
                    << " " << _kf_max2
                    << " " << _kf_min3 << endl;
  
  // ---- Define maximum kf above which we neglect data.
  //      At times there are a lot of very high frequency kmers that
  //      interfere with the genome size estimate. Those kmers 
  //      need to be thrown away.
  //
  //      Just a rough estimate for kf_hi is good enough.
  //      Assume:
  //  
  //          ndk[kf >= 2 . cov]  ~=  4 * ndk[2 . cov] / kf^2
  //         
  //      then lets find kf_hi for which:
  // 
  //          ndk[kf_hi] = 1   =>  kf_hi = sqrt(4 * ndk[2 . cov])
  //
  //
    
    
  _kf_hi = (2 * _kf_max2 < ndk.size()) ? 
    _kf_max2 * sqrt(4 * ndk[2 * _kf_max2] * _kf_max2) :
    _kf_max2 * sqrt(4 * ndk.back()        * _kf_max2);      
      
  if (_kf_hi > kf_ceil) 
    _kf_hi = kf_ceil;
    
  if (verbosity > 1) cout << "kf_hi           = " << _kf_hi << endl;

  // ---- number of read kmers in the various intervals

  _nk_total       = cnk.back();
  _nk_bad_low_kf  = cnk[_kf_min1];
  _nk_good_snp    = cnk[_kf_min2] - cnk[_kf_min1];
  _nk_good_uniq   = cnk[_kf_min3] - cnk[_kf_min2];
  _nk_good_rep    = cnk[_kf_hi]   - cnk[_kf_min3];
  _nk_bad_high_kf = _nk_total     - cnk[_kf_hi];
     
  if (verbosity > 1) cout << "nk_total        = " << _nk_total       << endl;
  if (verbosity > 1) cout << "nk_bad_low_kf   = " << _nk_bad_low_kf  << endl;
  if (verbosity > 1) cout << "nk_good_snp     = " << _nk_good_snp    << endl;
  if (verbosity > 1) cout << "nk_good_uniq    = " << _nk_good_uniq   << endl;
  if (verbosity > 1) cout << "nk_good_rep     = " << _nk_good_rep    << endl;
  if (verbosity > 1) cout << "nk_bad_high_kf  = " << _nk_bad_high_kf << endl;

    
  // ---- number of distinct kmers in the various intervals

  _ndk_total       = cndk.back();
  _ndk_bad_low_kf  = cndk[_kf_min1];
  _ndk_good_snp    = cndk[_kf_min2] - cndk[_kf_min1];
  _ndk_good_uniq   = cndk[_kf_min3] - cndk[_kf_min2];
  _ndk_good_rep    = cndk[_kf_hi]   - cndk[_kf_min3];
  _ndk_bad_high_kf = _ndk_total     - cndk[_kf_hi];
     
  if (verbosity > 1) cout << "ndk_total       = " << _ndk_total       << endl;
  if (verbosity > 1) cout << "ndk_bad_low_kf  = " << _ndk_bad_low_kf  << endl;
  if (verbosity > 1) cout << "ndk_good_snp    = " << _ndk_good_snp    << endl;
  if (verbosity > 1) cout << "ndk_good_uniq   = " << _ndk_good_uniq   << endl;
  if (verbosity > 1) cout << "ndk_good_rep    = " << _ndk_good_rep    << endl;
  if (verbosity > 1) cout << "ndk_bad_high_kf = " << _ndk_bad_high_kf << endl;


  // ---- kmer coverage C_K
 
  _kf_ave_uniq = double(_nk_good_uniq) / double(_ndk_good_uniq);
  if (verbosity > 1) cout << "kf_ave_uniq     = " << _kf_ave_uniq << endl;

  // ---- base coverage C_1
  // C_1 = C_K * L / ( L - K + 1)
  
  _coverage = (read_len > _K) ? _kf_ave_uniq * read_len / (read_len - _K + 1) : 0;
  if (verbosity > 1) cout << "read_len        = " << read_len << endl;
  if (verbosity > 1) cout << "coverage        = " << _coverage << endl;

  // ---- genome size
  _genome_size            = (_nk_total - _nk_bad_low_kf - _nk_bad_high_kf) / _kf_ave_uniq;

  _genome_size_unique     = _ndk_good_uniq + _ndk_good_snp / 2;


  if (false) {
    _genome_size_repetitive = size_t(0.5 + double(_nk_good_rep) / _kf_ave_uniq);

    _genome_size            = _genome_size_unique + _genome_size_repetitive;
  }
  else {
    _genome_size_repetitive = _genome_size - _genome_size_unique;
  }

  // ---- validate genome size estimate

  if (_nk_total / _genome_size < _coverage_min) {
    _genome_size = 0;
    _genome_size_unique = 0;
    _genome_size_repetitive = 0;
    _d_SNP = 0;
    _stddev_bias = 0;
  }
  else {

    // ---- SNP rate estimation
    //      Now, here's a funky expression for inter-SNP distance:
    //     
    //         d_SNP = 1 / ( 1 - (1 - ndk_good_snp / 2 G)^(1/K) )
    //
    //      It assumes uniform distribution of SNPs over the genome.
    //      It accounts for the reduction in SNP kmer counts when polymorphism rate is very high
    //      because of SNPs that are closer than K. 
    //      For low polymorphism rate the expression reduces to 
    //
    //         d_SNP = 2 K G / ndk_good_snp = G / n_snps
    //
    //      The factor of 2 accounts for the factor that ndk_good_snp is the count of 
    //      all kmers with SNPs so it equals: 2 K n_snps
    if (ploidy == 2) {
      _d_SNP = (_ndk_good_snp > 0) ?
        1.0 / (1.0 - pow(1.0 - 0.5 * _ndk_good_snp / _genome_size, 1.0 / _K)) : 1000000;
      
      if (verbosity > 1) cout << "SNP distance    = " << _d_SNP << endl;
    }


    // ---- Bias estimation
    {

      double sig2_SNP = 0;
      if (ploidy == 2) {
        // first, compute the variance of a Bernoulli distribution associated with SNPs and CN=1
        // with average at 1 and defined at x as the SNPs and at 2x as the CN=1 kmers
        // work out the formulas on your own to convince yourself that this is correct
        //
        double alpha = double(_ndk_good_snp) / double(_ndk_good_uniq);
        double tmp = 1.0 / (alpha + 2.0);
        sig2_SNP = alpha * tmp * tmp; 
        if (verbosity > 1) cout << "alpha           = " << setprecision(5) << alpha << endl;
        if (verbosity > 1) cout << "sig2 SNP        = " << setprecision(5) << sig2_SNP << endl;
        if (verbosity > 1) cout << "sig SNP         = " << setprecision(5) << sqrt(sig2_SNP) << endl;
      }

      
      // second, compute the sample mean and sample variance of the data
      // between min1 and min3
      double sum = 0;
      double mu = 0;
      double sig2 = 0;
      double sig = 0;
      while (sig == 0 || _kf_min3 - mu < 4 * sig || _kf_min3 < 3 * mu) {
	if (sig > 0)
	  _kf_min3++;

	sum = 0;
	for (size_t kf = _kf_min1; kf != _kf_min3; kf++) sum += ndk[kf];
	mu = 0;
	for (size_t kf = _kf_min1; kf != _kf_min3; kf++) mu += kf * ndk[kf];
	mu /= sum;

	sig2 = 0;
	for (size_t kf = _kf_min1; kf != _kf_min3; kf++) sig2 += (kf - mu) * (kf - mu) * ndk[kf];
	sig2 /= sum;
      
	sig = sqrt(fabs(sig2));
	if (verbosity > 1) cout << "kf_min3         = " << _kf_min3 << endl;
	if (verbosity > 1) cout << "mu              = " << mu << endl;
	if (verbosity > 1) cout << "sig             = " << sig << endl;
      }
            
      if (verbosity > 1) cout << "kf mean          = " << setprecision(5) << mu << endl;
      if (verbosity > 1) cout << "kf sig2          = " << setprecision(5) << sig2 << endl;
      if (verbosity > 1) cout << "kf sig           = " << setprecision(5) << sig << endl;
      if (verbosity > 1) cout << "kf sig/mean      = " << setprecision(5) << sig / mu << endl;

      // finally, compute the bias variance not forgetting to subtract 
      // the variance associated with SNPs
      double sig2_bias = (sig2 - mu)/(mu * mu) - sig2_SNP;
      _stddev_bias = sqrt(sig2_bias);
      if (verbosity > 1) cout << "sig2_bias       = " << setprecision(5) << sig2_bias << endl;
      if (verbosity > 1) cout << "stddev_bias     = " << setprecision(5) << _stddev_bias << endl;

    }
  }
}






void KmerSpectrum::to_text_file(const String & head,
                                const unsigned verbosity) const 
{
  const size_t nkf = this->size();

  size_t ndk_total = 0; // n distinct kmers
  size_t nk_total = 0;  // n kmers
  for (size_t kf = 0; kf != nkf; kf++) {
    ndk_total += (*this)[kf];
    nk_total  += (*this)[kf] * kf;
  }
    
  // ---- writing kmer frequency stats to disk
  const String fn = head_K(head) + ".kspec";
  if (verbosity) 
    cout << tag() << "Writing kmer spectrum to '" << fn << "'." << endl;
    
  ofstream outfs;
  outfs.open(fn.c_str());
  outfs << ("#"
            " 1:kmer_frequency"
            " 2:num_distinct_kmers"
            " 3:2/total_num_distinct_kmers"
            " 4:cummulative(2)"
            " 5:4/total_num_distinct_kmers"
            " 6:1*2"
            " 7:6/total_num_kmers"
            " 8:cummulative(6)"
            " 9:8/total_num_kmers") 
        << endl;

  size_t cndk = 0;
  size_t cnk = 0;
  for (size_t kf = 0; kf != nkf; kf++) {
    const size_t ndk = (*this)[kf];
    if (ndk > 0) {
      const size_t nk  = ndk * kf;
	
      cndk += ndk;
      cnk  += nk;
	
      outfs << " " << kf
            << " " << ndk
            << " " << setprecision(10) << double(ndk) / double(ndk_total) 
            << " " << cndk 
            << " " << setprecision(10) << double(cndk) / double(ndk_total) 
            << " " << nk
            << " " << setprecision(10) << double(nk) / double(nk_total) 
            << " " << cnk 
            << " " << setprecision(10) << double(cnk) / double(nk_total) 
            << endl;
    }
  }
  outfs.close();
}


void KmerSpectrum::from_text_file(const String & head,
                                  const unsigned verbosity)
{
  // ---- reading kmer frequency stats from disk
  const String fn = head_K(head) + ".kspec";
  if (verbosity) 
    cout << tag() << "Reading kmer spectrum from '" << fn << "'." << endl;
    
  this->clear();
  char buf[256];

  ifstream infs;
  infs.open(fn.c_str());
  size_t kf_prev;
  size_t kf = 0;
  size_t ndk;
  while (infs.good()) {
    infs.getline(buf, 256);
    if (kf == 0) {
      kf++;
      ndk = 0;
    }
    else {
      kf_prev = kf;
      sscanf(buf, "%lu %lu", &kf, &ndk);
    }
    this->resize(kf + 1, 0);
    (*this)[kf] = ndk;
  }
  infs.close();
}





















void KmerQualitySpectra::to_text_file(const String & head) const 
{
  const size_t nq = (*this).size();
  const size_t nkf = (*this)[0].size();
    

  // ---- find the total number of kmers at each quality score
  vec<size_t> ndk_total(nq, 0);
  vec<size_t> nk_total(nq, 0);
  for (size_t q = 0; q != nq; q++) {
    for (size_t kf = 0; kf != nkf; kf++) {
      ndk_total[q] += (*this)[q][kf];
      nk_total[q]  += (*this)[q][kf] * kf;
    }
  }


  // ---- open file
  const String fn = head + "." + ToString(_K) + "mer.kqspec";
  ofstream os;
  os.open(fn.c_str());


  {
    os << "# index [0]: 1:quality  2-5:num_distinct_kmers 6-9:num_kmers" << endl;
    size_t cndk = 0;
    size_t cnk = 0;
    for (size_t q = 1; q != nq; q++) {
      cndk += ndk_total[q];
      cnk  += nk_total[q];
      os << setw(12) << q
         << " " << setw(12) << ndk_total[q]
         << " " << setw(12) << double(ndk_total[q]) / double(ndk_total[0])
         << " " << setw(12) << cndk
         << " " << setw(12) << double(cndk) / double(ndk_total[0])
         << " " << setw(12) << nk_total[q]
         << " " << setw(12) << double(nk_total[q]) / double(nk_total[0])
         << " " << setw(12) << cnk
         << " " << setw(12) << double(cnk) / double(nk_total[0])
         << endl;
    }
    os << endl << endl;
  }

  {
    os << "# index [1]: 1:kmer_freq 2-...:num_distinct_kmers[q]" << endl;
    for (size_t kf = 0; kf != nkf; kf++) {
      if ((*this)[0][kf] != 0) {
        os << setw(12) << kf;

        for (size_t q = 0; q != nq; q++)  
          os << " " << setw(12) << (*this)[q][kf];

        os << endl;
      }
    }
    os << endl << endl;
  }


  {
    os << "# index [2]: 1:kmer_freq 2-...:num_distinct_kmers[q] / total_num_distinct_kmers" << endl;
    for (size_t kf = 0; kf != nkf; kf++) {
      if ((*this)[0][kf] != 0) {
        os << setw(12) << kf;

        for (size_t q = 0; q != nq; q++) 
          os << " " << setw(12) 
             << ((ndk_total[q]) ? double((*this)[q][kf]) / double(ndk_total[q]) : 0);

        os << endl;
      }
    }
    os << endl << endl;
  }


  {
    os << "# index [3]: 1:kmer_freq 2-...:cumulative_num_distinct_kmers[q]" << endl;
    vec<size_t> cndk(nq, 0);
    for (size_t kf = 0; kf != nkf; kf++) {
      if ((*this)[0][kf] != 0) {
        os << setw(12) << kf;
	  
        for (size_t q = 0; q != nq; q++) {
          cndk[q] += (*this)[q][kf];
          os << " " << setw(12) << cndk[q];
        }

        os << endl;
      }
    }
    os << endl << endl;
  }
    

  {
    os << "# index [4]: 1:kmer_freq 2-...:cumulative_num_distinct_kmers[q] / total_num_distinct_kmers" << endl;
    vec<size_t> cndk(nq, 0);
    for (size_t kf = 0; kf != nkf; kf++) {
      if ((*this)[0][kf] != 0) {
        os << setw(12) << kf;
	  
        for (size_t q = 0; q != nq; q++) {
          cndk[q] += (*this)[q][kf];
          os << " " << setw(12)
             << ((ndk_total[q]) ? double(cndk[q]) / double(ndk_total[q]) : 0);
        }

        os << endl;
      }
    }
    os << endl << endl;
  }



  {
    os << "# index [5]: 1:kmer_freq 2-...:num_kmers[q] = kmer_freq * num_distinct_kmers[q]" << endl;
    for (size_t kf = 0; kf != nkf; kf++) {
      if ((*this)[0][kf] != 0) {
        os << setw(12) << kf;

        for (size_t q = 0; q != nq; q++)  
          os << " " << setw(12) << kf * (*this)[q][kf];

        os << endl;
      }
    }
    os << endl << endl;
  }


  {
    os << "# index [6]: 1:kmer_freq 2-...:num_kmers[q] / total_num_kmers" << endl;
    for (size_t kf = 0; kf != nkf; kf++) {
      if ((*this)[0][kf] != 0) {
        os << setw(12) << kf;

        for (size_t q = 0; q != nq; q++)  
          os << " " << setw(12) 
             << ((nk_total[q]) ? double(kf * (*this)[q][kf]) / double(nk_total[q]) : 0);

        os << endl;
      }
    }
    os << endl << endl;
  }


  {
    os << "# index [7]: 1:kmer_freq 2-...:cumulative_num_kmers[q]" << endl;
    vec<size_t> cnk(nq, 0);
    for (size_t kf = 0; kf != nkf; kf++) {
      if ((*this)[0][kf] != 0) {
        os << setw(12) << kf;
	  
        for (size_t q = 0; q != nq; q++) {
          cnk[q] += kf * (*this)[q][kf];
          os << " " << setw(12) << cnk[q];
        }

        os << endl;
      }
    }
    os << endl << endl;
  }


    

  {
    os << "# index [8]: 1:kmer_freq 2-...:cumulative_num_kmers[q] / total_num_kmers" << endl;
    vec<size_t> cnk(nq, 0);
    for (size_t kf = 0; kf != nkf; kf++) {
      if ((*this)[0][kf] != 0) {
        os << setw(12) << kf;
	  
        for (size_t q = 0; q != nq; q++) {
          cnk[q] += kf * (*this)[q][kf];
          os << " " << setw(12) 
             << ((nk_total[q]) ? double(cnk[q]) / double(nk_total[q]) : 0);
        }

        os << endl;
      }
    }
    os << endl << endl;
  }


  os.close();


}

