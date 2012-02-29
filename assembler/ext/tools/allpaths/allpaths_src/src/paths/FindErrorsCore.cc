///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

//  Refactored by: Filipe Ribeiro - Nov 2011 - <crdhelp@broadinstitute.org>

#include "Basevector.h"
#include "Qualvector.h"
#include "feudal/QualNibbleVec.h"
#include "Bitvector.h"
#include "FeudalMimic.h"

#include "reporting/PerfStat.h"

#include "paths/FindErrorsCore.h"

#include "kmers/KmerSpectra.h"

#include "kmers/naif_kmer/NaifKmerizer.h"

#include "kmers/naif_kmer/KernelErrorFinder.h"
#include "kmers/naif_kmer/KernelPreCorrector.h"


static inline 
String Tag(String S = "FEC") { return Date() + " (" + S + "): "; } 






// ---- PreCorrect algorithm
//
// Corrects obvious substitution errors in reads.
// The Precorrect module examines all 25-mers in all reads.  It performs a
// large, external sort in such a way as to bring together for consideration
// all the 25-mers that agree about the initial and final 12-mers, ignoring the
// central base.  For the each pile of 25-mers that agree about the flanking
// bases, and potentially disagree only about that central base, it examines
// the quality scores associated with the base call of the central base.  Piles
// of fewer than 6 are ignored.  For larger piles, the quality scores are
// summed separately for each of the four potential calls of the central base.
// We declare a winner among the four calls as the call having the greatest sum
// of quality scores, but only if the sum is 60 or greater.  If we have a
// winner, then we can also check for losers.  A loser call has no more than
// one call of quality 20 or more, and its sum of quality scores must be less
// than 1/4 that of the winner.  Reads containing loser-calls are candidates for
// having those calls corrected to the winning call.  One further round of
// sorting and screening applies the correction to a read only if it is well
// isolated from other proposed corrections:  no changes are made to either
// call when two proposed corrections are within 12 bases of each other (i.e.,
// when the flanking bases that built a pile are themselves suspect).  Whenever
// a correction is adopted, the associated quality score for that call is set
// to 0.
  

template <class QVV_t>
void pre_correct_parallel(const PC_Params  & pcp, 
                          const unsigned     K, 
                          BaseVecVec       * bases_p, 
                          QVV_t            * quals_p, 
                          KmerSpectrum     * kspec_p, 
                          const unsigned     VERBOSITY,
                          const unsigned     NUM_THREADS, 
                          const size_t       mem_mean_ceil) 
{
  typedef typename QVV_t::value_type QV_t;

  cout << Tag() << "Duplicating 'bases' into 'bases_new'." << endl;
  
  BaseVecVec bases_new = *bases_p;

  cout << Tag() << "Finding recommendations." << endl;

  ForceAssertLe(K, 29u);
  if (K <= 29) {
    PreCorrector<Kmer29H, QVV_t> pre_corrector(pcp, K, *bases_p, *quals_p, 
                                               & bases_new, kspec_p, NUM_THREADS);
    naif_kmerize(&pre_corrector, NUM_THREADS, VERBOSITY, mem_mean_ceil);
  }

  // ---- Make corrections (unless they are too close to each other)
   
  cout << Tag() << "Making corrections." << endl;

  const size_t n_reads = bases_p->size();
  size_t n_bases = 0;
  size_t n_reads_corr = 0;
  size_t n_corr_skipped = 0;
  size_t n_corr_total = 0;
  for (size_t ibv = 0; ibv < n_reads; ibv++) {
    BaseVec & bv     = (*bases_p)[ibv];
    QV_t    & qv     = (*quals_p)[ibv];
    BaseVec & bv_new = bases_new[ibv];
    const size_t nb = bv.size();
    n_bases += nb;
    
    // -- print out requested reads
    
    if (pcp.i_reads_print.count(ibv) > 0) {
      cout << setw(12) << ibv << endl;
      cout << hieroglyphs(bv) << endl;
      for (size_t ib = 0; ib < nb; ib++) 
        cout << ((bv[ib] != bv_new[ib]) ? hieroglyph(bv_new[ib]) : ToString(" "));
      cout << endl;     
    }


    // -- store all indexes of corrections

    vec<size_t> ibs_corr;
    for (size_t ib = 0; ib < nb; ib++) 
      if (bv[ib] != bv_new[ib])
        ibs_corr.push_back(ib);
    const size_t n_corr = ibs_corr.size();
      
    // -- skip the ones that are closer than K/2

    vec<bool> skip(n_corr, false);
    for (size_t i = 1; i < n_corr; i++)
      if (ibs_corr[i] - ibs_corr[i-1] <= K/2) 
        skip[i] = skip[i-1] = true;

    // -- apply valid corrections 

    size_t n_corr_applied = 0;
    for (size_t i = 0; i < n_corr; i++)
      if (!skip[i]) {
        size_t ib = ibs_corr[i];
        bv.set(ib, bv_new[ib]);
	qv.set(ib, 0);
        n_corr_applied++;
      }

    if (n_corr_applied > 0) {
      n_reads_corr++;
      n_corr_total += n_corr_applied;
    }
    n_corr_skipped += n_corr - n_corr_applied;

    // -- print out requested reads

    if (pcp.i_reads_print.count(ibv) > 0) 
      cout << hieroglyphs(bv) << endl << endl;
  }


  // ---- Write correction statistics

  cout << Tag() << setw(12) << n_reads << " reads." << endl; 
  cout << Tag() << setw(12) << n_reads_corr << " reads corrected (" 
       << fixed << setprecision(1) << 100.0 * n_reads_corr / n_reads << " %)." << endl;

  cout << Tag() << setw(12) << n_corr_total << " total corrections." << endl;
  cout << Tag() << setw(12) << fixed << setprecision(1) << double(n_corr_total) / double(n_reads_corr) 
       << " corrections per corrected read." << endl;
  cout << Tag() << setw(12) << n_corr_skipped << " corrections skipped." << endl;


  PerfStat::log() << std::fixed << std::setprecision(1)
                  << PerfStat("frac_bases_pre_corrected", 
                              "% of bases pre-corrected", 
                              100.0 * float(n_corr_total) / float(n_bases));
}



// Instantiate QualVec and QualNibbleVec

template void pre_correct_parallel(const PC_Params  & pcp, 
                                   const unsigned     K_PC, 
                                   BaseVecVec       * bases_p, 
                                   QualVecVec       * quals_p, 
                                   KmerSpectrum     * kspec_p, 
                                   const unsigned     VERBOSITY,
                                   const unsigned     NUM_THREADS, 
                                   const size_t       mem_mean_ceil);

template void pre_correct_parallel(const PC_Params  & pcp, 
                                   const unsigned     K_PC, 
                                   BaseVecVec       * bases_p, 
                                   QualNibbleVecVec * quals_p, 
                                   KmerSpectrum     * kspec_p, 
                                   const unsigned     VERBOSITY,
                                   const unsigned     NUM_THREADS, 
                                   const size_t       mem_mean_ceil);
















// ------ Capping quality scores in parallel usiong worklist


template<class QVV_t>
class QualCapperProc
{
  QVV_t            * _p_quals;
  const size_t       _QUAL_CEIL_RADIUS;
  const size_t       _n_threads;

  typedef typename QVV_t::value_type QV_t;

public:
  QualCapperProc(QVV_t         * p_quals, 
                 const size_t    QUAL_CEIL_RADIUS,
                 const size_t    n_threads)
    : _p_quals(p_quals),
      _QUAL_CEIL_RADIUS(QUAL_CEIL_RADIUS),
      _n_threads(n_threads)
  {}   

  void operator()(const size_t i_thread) 
  {
    const size_t n_qv = _p_quals->size();
    const size_t i0_qv = ( i_thread      * n_qv) / _n_threads;
    const size_t i1_qv = ((i_thread + 1) * n_qv) / _n_threads;

    vec<size_t> qv_new;
    for (size_t i_qv = i0_qv; i_qv != i1_qv; i_qv++) {
      QV_t & qv = (*_p_quals)[i_qv];
      const size_t n_q = qv.size();
      if (qv_new.size() < n_q) qv_new.resize(n_q);
      
      for (size_t i_q = 0; i_q < n_q; i_q++) {
        size_t q_min = qv[i_q];
        const size_t j_q0 = (i_q <=           _QUAL_CEIL_RADIUS) ?       0 : i_q - _QUAL_CEIL_RADIUS;
        const size_t j_q1 = (i_q >= n_q - 1 - _QUAL_CEIL_RADIUS) ? n_q - 1 : i_q + _QUAL_CEIL_RADIUS;
        for (size_t j_q = j_q0; j_q <= j_q1; j_q++) 
          if (q_min > qv[j_q]) 
            q_min = qv[j_q];
        qv_new[i_q] = q_min;
      }
      for (size_t i_q = 0; i_q < n_q; i_q++)
        qv.set(i_q, qv_new[i_q]);
    }
  }
};



template <class QVV_t>
void cap_quality_scores(QVV_t        * p_quals,
			const size_t   QUAL_CEIL_RADIUS,
                        const size_t   NUM_THREADS)
			
{
  QualCapperProc<QVV_t> capper(p_quals, QUAL_CEIL_RADIUS, NUM_THREADS);
  Worklist<size_t, QualCapperProc<QVV_t> > worklist(capper, NUM_THREADS - 1);
  
  for (size_t i = 0; i != NUM_THREADS; i++) 
    worklist.add(i);
}





// ---- class of counters to manage correction statistics


class CounterMap : public map<size_t, size_t>
{
public:
  size_t total() const 
  {
    size_t total = 0;
    for (CounterMap::const_iterator it = (*this).begin(); 
         it != (*this).end(); it++)
      total += it->second;
    return total;
  }

  void to_text_file(const String & filename) const
  {
    double tot = static_cast<double>(total());

    ofstream ofs;
    ofs.open(filename.c_str());

    ofs << "# 1:x 2:y 3:y/total 4:cum(y) 5:cum(y)/total" << endl;

    size_t cum = 0;
    for (CounterMap::const_iterator it = (*this).begin(); 
       it != (*this).end(); it++) {
      cum += it->second;
      ofs << " " << it->first
          << " " << it->second
          << " " << setprecision(10) << static_cast<double>(it->second) / tot
          << " " << cum
          << " " << setprecision(10) << static_cast<double>(cum) / tot
          << endl;
    }
    ofs.close();
  }
  /*
  CounterMap & operator += (const CounterMap & m) 
  {
    for (CounterMap::const_iterator it = m.begin(); it != m.end(); it++)
      (*this)[(*it).first] += (*it).second;
    return *this;
  } 
  */
};




class EF_Stats
{
public:
  size_t n_bases;
  size_t n_confirmed;
  size_t n_corrections;
  size_t n_conflicts;
  CounterMap n_corrections_pos;
  CounterMap n_corrections_quals;
  CounterMap n_conflicts_pos;
  CounterMap n_conflicts_quals;
  
  EF_Stats() : n_bases(0), n_confirmed(0), n_corrections(0), n_conflicts(0) {}
};





// ---- Apply all the corrections and compute correction statistics

template<class QVV_t>
void apply_corrections(const BaseVecVec & bases_new, 
                       const BitVecVec  & base_locked, 
                       BaseVecVec       * bases_p, 
                       QVV_t            * quals_p,
                       EF_Stats         * ef_stats_p)
{
  const size_t n_reads = bases_p->size();

  ef_stats_p->n_bases = 0;
  ef_stats_p->n_corrections = 0;
  ef_stats_p->n_conflicts = 0;
  ef_stats_p->n_confirmed = 0;
    
  CounterMap n_corrections_pos;
  CounterMap n_corrections_quals;
  
  
  for (size_t i_read = 0; i_read != n_reads; i_read++) {
    
    const size_t n_bases = (*bases_p)[i_read].size();
    ef_stats_p->n_bases += n_bases;

    for (size_t i_base = 0; i_base < n_bases; i_base++) {

      if (bases_new[i_read][i_base] != (*bases_p)[i_read][i_base]) {  
        const unsigned q = (*quals_p)[i_read][i_base];

        if (! base_locked[i_read][i_base]) {

          // process only base locations that are NOT base_locked.

          (*bases_p)[i_read].set(i_base, bases_new[i_read][i_base]);
          (*quals_p)[i_read].set(i_base, 0);
          
          ef_stats_p->n_corrections++;
          ef_stats_p->n_corrections_pos[i_base]++;
          ef_stats_p->n_corrections_quals[q]++;
        }
        else {
          ef_stats_p->n_conflicts++;
          ef_stats_p->n_conflicts_pos[i_base]++;
          ef_stats_p->n_conflicts_quals[q]++;
        }
      }
      else {
        if (base_locked[i_read][i_base])
          ef_stats_p->n_confirmed++;
      }
    }
  }
 
  
}



// ---- Find correction recommendations in parallel

template<class QVV_t>
void find_errors_parallel(const EF_Params  & efp, 
                          const size_t       K, 
                          BaseVecVec       * bases_p, 
                          QVV_t            * quals_p, 
                          const size_t       NUM_CYCLES,
                          const unsigned     VERBOSITY,
                          const unsigned     NUM_THREADS, 
                          const size_t       mem_mean_ceil)
{

  for (size_t i_cycle = 0; i_cycle != NUM_CYCLES; i_cycle++) {
    if (VERBOSITY > 0) cout << Tag() << "cycle = " << i_cycle << "/" << NUM_CYCLES << endl;
    

    // ---- cap quality scores
    if (efp.qual_ceil_radius > 0 && (i_cycle == 0 || efp.qcr_always)) {
      if (VERBOSITY > 0) cout << Tag() << "Renormalizing quality scores." << endl;
      cap_quality_scores(quals_p, efp.qual_ceil_radius, NUM_THREADS);
    }



    // ---- Keep track of which bases have been recommended for corrections,
    //      and which bases have been base_locked.
    if (VERBOSITY > 0) cout << Tag() << "Duplicating 'bases' into 'bases_new'." << endl;
    BaseVecVec bases_new = *bases_p;
      
    if (VERBOSITY > 0) cout << Tag() << "Mimicking 'bases' into 'base_locked'." << endl;
    BitVecVec base_locked;
    Mimic(*bases_p, base_locked);

    // ---- PARALLEL PROCESSING
    //      Each parcel is processed at a time and the correction recommendations 
    //      are stored in 'bases_new' and 'base_locked'.
      
    if (VERBOSITY > 0) cout << Tag() << "Gathering list of base-locations to correct." << endl;
      
    ForceAssertLe(K, 124u);
    if (K <= 29) {
      ErrorFinder<Kmer29, QVV_t> error_finder(efp, K, *bases_p, *quals_p, 
                                              &bases_new, &base_locked, NUM_THREADS);
      naif_kmerize(&error_finder, NUM_THREADS, VERBOSITY, mem_mean_ceil);
    }
    else if (K <= 60) {
      ErrorFinder<Kmer60, QVV_t> error_finder(efp, K, *bases_p, *quals_p, 
                                              &bases_new, &base_locked, NUM_THREADS);
      naif_kmerize(&error_finder, NUM_THREADS, VERBOSITY, mem_mean_ceil);
    }
    else if (K <= 124) {
      ErrorFinder<Kmer124, QVV_t> error_finder(efp, K, *bases_p, *quals_p, 
                                               &bases_new, &base_locked, NUM_THREADS);
      naif_kmerize(&error_finder, NUM_THREADS, VERBOSITY, mem_mean_ceil);
    }

      
    // ---- Apply corrections
    if (VERBOSITY > 0) cout << Tag() << "Applying the recommended corrections." << endl;

    EF_Stats ef_stats;
    apply_corrections(bases_new, base_locked, bases_p, quals_p, &ef_stats);


    // ---- Report results 
    if (VERBOSITY > 0) {
      cout << Tag() << "n_bases       = " << setw(12) << ef_stats.n_bases << endl;
      cout << Tag() << "n_confirmed   = " << setw(12) << ef_stats.n_confirmed << endl;
      cout << Tag() << "n_corrections = " << setw(12) << ef_stats.n_corrections << endl;
      cout << Tag() << "n_conflicts   = " << setw(12) << ef_stats.n_conflicts << endl;
      
      PerfStat::log() << std::fixed << std::setprecision(1)
                      << PerfStat("frac_bases_confirmed_" + ToString(i_cycle), 
                                  "% of bases confirmed in cycle " + ToString(i_cycle), 
                                  100.0 * float(ef_stats.n_confirmed) / float(ef_stats.n_bases));
      PerfStat::log() << std::fixed << std::setprecision(2)
                      << PerfStat("frac_bases_corrected_" + ToString(i_cycle), 
                                  "% of bases corrected in cycle " + ToString(i_cycle), 
                                  100.0 * float(ef_stats.n_corrections) / float(ef_stats.n_bases));
      PerfStat::log() << std::fixed << std::setprecision(2)
                      << PerfStat("frac_conflicts_" + ToString(i_cycle), 
                                  "% of bases with conflicting corrections in cycle " + ToString(i_cycle), 
                                  100.0 * float(ef_stats.n_conflicts) / float(ef_stats.n_bases));

      /* 
        const String head_c = out_full_edit_head + "." + ToString(K) + "mer.cycle" + ToString(i_cycle);
        ef_stats.n_corrections_pos.to_text_file(head_c + ".n_corrs_pos.count.dat");
        ef_stats.n_corrections_quals.to_text_file(head_c + ".n_corrs_quals.count.dat");
        ef_stats.n_conflicts_pos.to_text_file(head_c + ".n_confls_pos.count.dat");
        ef_stats.n_conflicts_quals.to_text_file(head_c + ".n_confls_quals.count.dat");
      */
    }
  
  } 


}

// ---- Instantiations

template void find_errors_parallel(const EF_Params  & efp, 
                                   const size_t       K, 
                                   BaseVecVec       * bases_p, 
                                   QualVecVec       * quals_p, 
                                   const size_t       NUM_CYCLES,
                                   const unsigned     VERBOSITY,
                                   const unsigned     NUM_THREADS, 
                                   const size_t       mem_mean_ceil);

template void find_errors_parallel(const EF_Params  & efp, 
                                   const size_t       K, 
                                   BaseVecVec       * bases_p, 
                                   QualNibbleVecVec * quals_p, 
                                   const size_t       NUM_CYCLES,
                                   const unsigned     VERBOSITY,
                                   const unsigned     NUM_THREADS, 
                                   const size_t       mem_mean_ceil);




// ---- Single thread call to find errors
//      Useful when correcting small local read sets in modules other than FindErrors


template<class QVV_t>
void find_errors(const EF_Params  & efp, 
                 const size_t       K, 
                 BaseVecVec       * bases_p, 
                 QVV_t            * quals_p, 
                 const size_t       NUM_CYCLES,
                 const size_t       VERBOSITY)
{
  const size_t n_threads = 1;
  const size_t mem_mean_ceil = 0;
  find_errors_parallel(efp, K, bases_p, quals_p, NUM_CYCLES, 
                       VERBOSITY, n_threads, mem_mean_ceil);
}




// Instantiate for QualVec and QualNibbleVec

template void find_errors(const EF_Params  & efp, 
                          const size_t       K, 
                          BaseVecVec       * bases_p, 
                          QualVecVec       * quals_p, 
                          const size_t       NUM_CYCLES,
                          const size_t       VERBOSITY);


template void find_errors(const EF_Params  & efp, 
                          const size_t       K, 
                          BaseVecVec       * bases_p, 
                          QualNibbleVecVec * quals_p, 
                          const size_t       NUM_CYCLES,
                          const size_t       VERBOSITY);




