///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
  ParcelKmersTools

  2009-07 Josh Burton (initial developer)
  2010-03 Filipe Ribeiro

  A toolkit of useful functions that use the kmers parcels paradigm.
 
*/

#ifndef KMER_PARCELS_TOOLS_H
#define KMER_PARCELS_TOOLS_H


#include "kmers/KmerParcels.h"



inline 
String KPTTag(String S = "PK") { return Date() + " (" + S + "): "; } 



inline bool GetNextMatchingKmerBatches(KmerParcelReader & parcel1_reader,
                                KmerParcelReader & parcel2_reader)
{
  // always read the next k1
  if (!parcel1_reader.GetNextKmerBatch()) return false;

  // always read the next k2
  if (!parcel2_reader.GetNextKmerBatch()) return false;


  const NaifBuffer * kmer1 = & parcel1_reader.CurrentKmerBatch().Kmer();
  const NaifBuffer * kmer2 = & parcel2_reader.CurrentKmerBatch().Kmer();

  while (1) {

    if (*kmer1 < *kmer2) { // parcel 1 must catch up 

      if (!parcel1_reader.GetNextKmerBatch()) return false;
      kmer1 = & parcel1_reader.CurrentKmerBatch().Kmer();
    }
    else if (*kmer1 > *kmer2) { // parcel 2 must catch up 
      
      if (!parcel2_reader.GetNextKmerBatch()) return false;
      kmer2 = & parcel2_reader.CurrentKmerBatch().Kmer();

    }
    else { // kmer1 == kmer2
      
      return true;
      
    }
    
  }
  
}
  




inline bool GetExclusiveKmerBatch(KmerParcelReader & parcel1_reader,
                           KmerParcelReader & parcel2_reader)
{
  size_t nb1 = parcel1_reader.GetNumKmerBatches();
  size_t nb2 = parcel2_reader.GetNumKmerBatches();

  // always read the next 1
  if (!parcel1_reader.GetNextKmerBatch()) return false;

  // if never read 2 read the first
  if (parcel2_reader.GetNumKmerBatchesRead() == 0)
    if (!parcel2_reader.GetNextKmerBatch()) return true;

  const NaifBuffer * kmer1 = & parcel1_reader.CurrentKmerBatch().Kmer();
  const NaifBuffer * kmer2 = & parcel2_reader.CurrentKmerBatch().Kmer();

  while (1) {
    //cout << "  -- while(1) ------" << endl;

    if (*kmer1 < *kmer2) { // => all done; kmer1 is exclusive

      return true;

    }
    else if (*kmer1 > *kmer2) { // => k2 must catch up

      if (!parcel2_reader.GetNextKmerBatch()) return true;

      kmer2 = & parcel2_reader.CurrentKmerBatch().Kmer();

    }
    else { // (kmer1 == kmer2) => skip both

      if (!parcel1_reader.GetNextKmerBatch()) return false;
      if (!parcel2_reader.GetNextKmerBatch()) return true;

      kmer1 = & parcel1_reader.CurrentKmerBatch().Kmer();
      kmer2 = & parcel2_reader.CurrentKmerBatch().Kmer();

    } 

  }
}






















// ----------------------------------------------
// CompareKmerParcels
// ----------------------------------------------



inline void CountMatchingKmersInParcels(KmerParcelReader & parcel1_reader,
                                 KmerParcelReader & parcel2_reader,
                                 ComparativeKmerFrequencyCounters * kmer_freqs)
{
  size_t K = parcel1_reader.GetK();
  NaifBuffer kmer1(2*K);
  NaifBuffer kmer2(2*K);
  size_t freq1, freq2;
  bool found1, found2;
  
  NaifTimer timer_all;

  
  timer_all.Start();
  
  found1 = parcel1_reader.GetNextKmer(&kmer1, &freq1);
  found2 = parcel2_reader.GetNextKmer(&kmer2, &freq2);
  
  while (found1 && found2) {
    
    if (kmer1 < kmer2) {  // parcel 1 has to catch up to parcel 2

      kmer_freqs->only1[freq1]++;
      found1 = parcel1_reader.GetNextKmer(&kmer1, &freq1);
    }
    else if (kmer1 > kmer2) {  // parcel 2 has to catch up to parcel 1

      kmer_freqs->only2[freq2]++;
      found2 = parcel2_reader.GetNextKmer(&kmer2, &freq2);
    }
    else {  // kmer1 and kmer2 match
      kmer_freqs->both1[freq1]++;
      found1 = parcel1_reader.GetNextKmer(&kmer1, &freq1);

      kmer_freqs->both2[freq2]++;
      found2 = parcel2_reader.GetNextKmer(&kmer2, &freq2);
    }
  }

  while (found1) {
    kmer_freqs->only1[freq1]++;
    found1 = parcel1_reader.GetNextKmer(&kmer1, &freq1);
  }

  while (found2) {
    kmer_freqs->only2[freq2]++;
    found2 = parcel2_reader.GetNextKmer(&kmer2, &freq2);
  }

  timer_all.Stop();

  double timer_read = parcel1_reader.GetTimer() + parcel2_reader.GetTimer();
  double timer_else = timer_all.GetTimer() - timer_read;

  /*
  if (_verbose)
    cout << (KPTTag() + "timers(sec):" + 
	     " read= "  + ToString(timer_read, 1)  + " " + ToString(100*timer_read/timer_all, 1)  + "%" + 
	     " else= "  + ToString(timer_else, 1)  + " " + ToString(100*timer_else/timer_all, 1)  + "%")
	 << endl;
  */
}




class CountMatchingKmersProc 
{
private:
  KmerParcelsStore & _parcels1;
  KmerParcelsStore & _parcels2;

  vec<ComparativeKmerFrequencyCounters> * _kmer_freqs_vec;

  NaiveThreadPool & _thread_pool;


public:
  CountMatchingKmersProc(KmerParcelsStore & parcels1,
                         KmerParcelsStore & parcels2,
                         vec<ComparativeKmerFrequencyCounters> * kmer_freqs_vec,
                         NaiveThreadPool & thread_pool)
    : 
      _parcels1(parcels1),
      _parcels2(parcels2),
      _kmer_freqs_vec(kmer_freqs_vec),
      _thread_pool(thread_pool)
  {}
    

  void operator() (size_t parcel_ID) {

    const size_t thread_ID = _thread_pool.AssignThread(parcel_ID);

    String tag = "CKPP";
    //if (_verbose)
      cout << KPTTag(tag) << "Processing parcel " << parcel_ID 
           << " of " << _parcels1.GetNumParcels()
           << " on thread " << thread_ID << endl;
    {
      KmerParcelReader parcel1_reader(_parcels1, parcel_ID);
      KmerParcelReader parcel2_reader(_parcels2, parcel_ID);
      

      CountMatchingKmersInParcels(parcel1_reader,
                                  parcel2_reader,
                                  &((*_kmer_freqs_vec)[thread_ID]));
    }

    _thread_pool.UnassignThread(parcel_ID);
  }

};









inline void CompareKmerParcels(KmerParcelsStore & parcels1,
                               KmerParcelsStore & parcels2,
                               ComparativeKmerFrequencyCounters * kmer_freqs, 
                               const size_t num_threads = 1)
{
  String tag = "CKPP";

  ForceAssertEq(parcels1.GetNumParcels(), parcels2.GetNumParcels());

  const size_t num_parcels = parcels1.GetNumParcels();

  vec<ComparativeKmerFrequencyCounters> kmer_freqs_vec(num_threads);

  {
    pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
    NaiveThreadPool thread_pool(num_threads, num_parcels, &lock);
    
    CountMatchingKmersProc proc(parcels1, parcels2, &kmer_freqs_vec, thread_pool);
    
  
    Worklist<size_t, CountMatchingKmersProc> 
      worklist(proc, num_threads - 1);
    
    for (size_t parcel_ID = 0; parcel_ID < num_parcels; parcel_ID++)
      worklist.add(parcel_ID);
      
  }

  typedef vec<ComparativeKmerFrequencyCounters>::const_iterator CKFCVecIter;
  for (CKFCVecIter it = kmer_freqs_vec.begin(); it != kmer_freqs_vec.end(); it++)
    (*kmer_freqs) += (*it);
      
}




// wrapper for function above that writes results to disk
inline void CompareKmerParcels(KmerParcelsDiskStore & parcels1,
                        KmerParcelsDiskStore & parcels2,
                        const String & label1,
                        const String & label2,
                        const size_t num_threads = 1) 
{
  String tag = "CKP";

  ComparativeKmerFrequencyCounters kmer_freqs;

  CompareKmerParcels(parcels1, parcels2, &kmer_freqs, num_threads);
                     
  cout << KPTTag(tag) << "Writing results" << endl;

  String Dir1Name = parcels1.GetDirectoryName();
  String Dir2Name = parcels2.GetDirectoryName();


  /*
  kmer_freqs.only1.Write(Dir1Name + "/" + label2 + ".not_present.kmer_frequencies.count.dat");
  kmer_freqs.both1.Write(Dir1Name + "/" + label2 + ".present.kmer_frequencies.count.dat");
  kmer_freqs.only2.Write(Dir2Name + "/" + label1 + ".not_present.kmer_frequencies.count.dat");
  kmer_freqs.both2.Write(Dir2Name + "/" + label1 + ".present.kmer_frequencies.count.dat");
  */

  size_t only1 = kmer_freqs.only1.Total();
  size_t only2 = kmer_freqs.only2.Total();
  size_t both1 = kmer_freqs.both1.Total();
  size_t both2 = kmer_freqs.both2.Total();
  size_t all1 = only1 + both1;
  size_t all2 = only2 + both2;
  

  cout << KPTTag(tag) << "GREPABLE: " 
       << setw(12) << all1
       << " (" << setw(5) << fixed  << setprecision(1) << 100.0 << "%)"
       << "  kmers in      '" << label1 << "'." << endl;
  cout << KPTTag(tag) << "GREPABLE: " 
       << setw(12) << both1
       << " (" << setw(5) << fixed  << setprecision(1) 
       << 100.0 * static_cast<double>(both1) / static_cast<double>(all1)
       << "%)"
       << "  kmers also in '" << label2 << "'." << endl;
  cout << KPTTag(tag) << "GREPABLE: " 
       << setw(12) << only1
       << " (" << setw(5) << fixed  << setprecision(1) 
       << 100.0 * static_cast<double>(only1) / static_cast<double>(all1)
       << "%)"
       << "  kmers not in  '" << label2 << "'." << endl;

  cout << endl;

  cout << KPTTag(tag) << "GREPABLE: " 
       << setw(12) << all2
       << " (" << setw(5) << fixed  << setprecision(1) << 100.0 << "%)"
       << "  kmers in      '" << label2 << "'." << endl;
  cout << KPTTag(tag) << "GREPABLE: " 
       << setw(12) << both2
       << " (" << setw(5) << fixed  << setprecision(1) 
       << 100.0 * static_cast<double>(both2) / static_cast<double>(all2)
       << "%)"
       << "  kmers also in '" << label1 << "'." << endl;
  cout << KPTTag(tag) << "GREPABLE: " 
       << setw(12) << only2
       << " (" << setw(5) << fixed  << setprecision(1) 
       << 100.0 * static_cast<double>(only2) / static_cast<double>(all2)
       << "%)"
       << "  kmers not in  '" << label1 << "'." << endl;
}












inline bool FindMatchingKmerBatch(const KmerParcelReader & parcel_reader,
                           const NaifBuffer & kmer0)
{
  // if never read any batch, read the first
  if (parcel_reader.GetNumKmerBatchesRead() == 0)
    if (!parcel_reader.GetNextKmerBatch()) return false;

  while (true) {
    const NaifBuffer * kmer = & parcel_reader.CurrentKmerBatch().Kmer(); 
    if      (kmer0  < *kmer)                    return false;
    else if (kmer0 == *kmer)                    return true;
    else if (!parcel_reader.GetNextKmerBatch()) return false;
  }
}






#endif
