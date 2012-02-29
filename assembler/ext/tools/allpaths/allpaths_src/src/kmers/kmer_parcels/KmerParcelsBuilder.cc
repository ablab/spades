///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#include "kmers/KmerParcels.h"
#include "system/WorklistN.h"


inline String Tag(String S = "KPs") { return Date() + " (" + S + "): "; } 



void PrintValueAndPercentage(ostream & out, const String & s, const size_t val, const double total)
{
  out << s 
      << setw(12) << val 
      << " " << setw(5) << fixed << setprecision(1) << (100.0 * val) / total << "%" 
      << endl;
}



/* key_digit_table: A static lookup table to convert two-base codes into octal
 * (base-8) digits, for use in generating parcel_IDs.
 *
 * A "two-base code" is a number in the range [0,16), calculated from two bases
 * at the beginning and end of a kmer.  It's like a two-digit number in base 4,
 * where A=0, C=1, G=2, T=3.  For example, if a kmer begins with G and and ends
 * with C, the two-base code is 4*(2) + (1) = 9.
 *
 * The two-base codes are mapped onto the octal digits so that the digit is
 * invariant under the RC operator.  For example, AG gives the same digit as
 * CT.  Hence a kmer will have the same parcel_ID however it is canonicalized.
 *
 * The key_digit_table is designed to ensure that a roughly equal number of
 * kmers are placed in each parcel, making the parcels similar in size.
 * This is slightly perturbed if n_parcels is not a power of 2.
 * It is more strongly affected if the input bases are GC-rich or AT-rich; this
 * particularly throws off digits '6' and '7'.  For example:
 * If GC-content is 50%, then '6' and '7' both come out of the key_digit_table
 * with frequency 12.5%, as do all the other digits.  But if GC-content is 60%,
 * then '6' appears with frequency 8%, and '7' appears with frequency 18%.
 *
 */

const char key_digit_table[] =
  { 6, 0, 1, 4,
    2, 7, 4, 1,
    3, 5, 7, 0,
    5, 3, 2, 6 };


// ---- Computes a number in [0,7] from the bases at 
//      positions 'base_ID' and 'base_ID' + 'K' - 1
//      of the BaseVec 'bv'. 
inline 
size_t ComputeOctalKeyFromBaseVec(const BaseVec & bv, 
                                  const size_t base_ID,
                                  const size_t K)
{
  return key_digit_table[(bv[base_ID] << 2) + bv[base_ID + K - 1]];
}


// ---- Conputes the parcel_ID of the 'K'-mer of BaseVec 'bv' 
//      at position 'base_ID'. 
inline 
size_t ComputeParcelIDFromBaseVec(const BaseVec & bv,
                                  const size_t base_ID,
                                  const size_t K,
                                  const size_t n_parcels)
{
  if (n_parcels <= 8) {
    ForceAssertGt(n_parcels, 0u);
    const size_t key = ComputeOctalKeyFromBaseVec(bv, base_ID, K);
    return (key * n_parcels) >> 3;   // ... / 8
  }
  else if (n_parcels <= 64) {
    const size_t key = ((ComputeOctalKeyFromBaseVec(bv, base_ID    , K    ) << 3) + 
                        (ComputeOctalKeyFromBaseVec(bv, base_ID + 1, K - 2)     ));
    return (key * n_parcels) >> 6;   // ... / 64
  }
  else if (n_parcels <= 512) {
    const size_t key = ((ComputeOctalKeyFromBaseVec(bv, base_ID    , K    ) << 6) + 
                        (ComputeOctalKeyFromBaseVec(bv, base_ID + 1, K - 2) << 3) +
                        (ComputeOctalKeyFromBaseVec(bv, base_ID + 2, K - 4)     ));
    return (key * n_parcels) >> 9;   // ... / 512
  }
  else if (n_parcels <= 4096) {
    const size_t key = ((ComputeOctalKeyFromBaseVec(bv, base_ID    , K    ) << 9) + 
                        (ComputeOctalKeyFromBaseVec(bv, base_ID + 1, K - 2) << 6) +
                        (ComputeOctalKeyFromBaseVec(bv, base_ID + 2, K - 4) << 3) +
                        (ComputeOctalKeyFromBaseVec(bv, base_ID + 3, K - 6)     ));
    return (key * n_parcels) >> 12;  // ... / 4096
  }
  else if (n_parcels <= 32768) {
    const size_t key = ((ComputeOctalKeyFromBaseVec(bv, base_ID    , K    ) << 12) + 
                        (ComputeOctalKeyFromBaseVec(bv, base_ID + 1, K - 2) <<  9) +
                        (ComputeOctalKeyFromBaseVec(bv, base_ID + 2, K - 4) <<  6) +
                        (ComputeOctalKeyFromBaseVec(bv, base_ID + 3, K - 6) <<  3) +
                        (ComputeOctalKeyFromBaseVec(bv, base_ID + 4, K - 8)      ));
    return (key * n_parcels) >> 15;  // ... / 32768
  }
  else {
    cout << "**** n_parcels = " << n_parcels << " > 32768!" << endl;
    cout << "     Expand ComputeParcelIDFromBaseVec(...) in kmers/kmer_parcels/KmerParcelsBuilder.cc"
      " to accomodate more parcels." << endl; 
  }
  ForceAssertLe(n_parcels, 32768ul);
  return 0;
}


// ---------------------------------------------
// KmerParcel<K> [Locked vec< KmerRecord<K> >]
// ---------------------------------------------


// several threads might add to this class independently while 
// processing a block of reads
template<size_t K>
void KmerParcel<K>::PushKmerRecords(const vec< KmerRecord<K> > & buf)
{
  Locker lock(*this);   // mutex protecting this block
  (*this).insert((*this).end(), 
                 buf.begin(), 
                 buf.end());
}
  


// ---- sort (*this) vector so that repeated kmers end up next to each other
//      ribeiro: Note: 
//      ribeiro:   first sorts over read_IDs and then over SIGNED pos
//      ribeiro:   meaning that, for the same read_ID, we get all the RCs kmers before the FWs
//      ribeiro:   DO WE NEED IT LIKE THIS? I DON'T KNOW/REMEMBER.  
template<size_t K>
void KmerParcel<K>::SortKmerParcel()
{
  _timer_sort.Start();
  sort((*this).begin(), (*this).end(), KmerRecord<K>::KmerReadIDSignedPosLT); 
  _timer_sort.Stop();
}



template<size_t K>
void KmerParcel<K>::WriteSortedKmerParcel(const size_t thread_ID, 
                                          const size_t parcel_ID, 
                                          KmerParcelsStore & parcels_store,
                                          KmerFrequencyStatistics & stats)
{
  String tag = "WSKP " + ToString(parcel_ID);

  // ---- collect all kmer locations, sort them, and write them to disk as a batch
  _timer_write.Start();
  
  // ---- get a writer for this parcel_ID
  KmerParcelWriter parcel_writer(parcels_store, parcel_ID);
  

  size_t n_batches = 0;

  vec<const KmerLoc *> kmer_locs_p;

  typedef typename vec< KmerRecord<K> >::const_iterator KmerRecordIter;
  for (KmerRecordIter it0 = (*this).begin(); 
       it0 != (*this).end(); it0 += kmer_locs_p.size()) {
    kmer_locs_p.clear();

    // ---- look for all the consecutive occurrences of this kmer
    KmerRecordIter it1 = it0;
    while (it1 != (*this).end() && KmerRecord<K>::KmersEQ(*it0, *it1)) {
      kmer_locs_p.push_back(static_cast<const KmerLoc*>(&(*it1)));
      ++it1;
    }
    const size_t kmer_freq = it1 - it0;

    // ---- update parcels on store (disk or memory)
    if (parcels_store.KmerValidator(kmer_freq))
      parcel_writer.PutKmerKmerLocs(it0->Kmer(), kmer_locs_p);

    // ---- update kmer freqeuncy stats
    stats.Update(kmer_locs_p, thread_ID, it0->GetKmerNumGC());
      
    n_batches++;
  }

  _timer_write.Stop();

  double t_sort = _timer_sort.GetTimer();
  double t_write = _timer_write.GetTimer();

  double t_all = t_sort + t_write;
  
  parcels_store.Log(Tag(tag) + 
		    "n_batches= "  + ToString(n_batches)  + 
		    " timers(sec): " + 
		    "sort= "  + ToString(t_sort, 1)  + 
		    " " + ToString(100*t_sort/t_all, 1)  + "% " + 
		    "write= " + ToString(t_write, 1) + 
		    " " + ToString(100*t_write/t_all, 1) + "%");
  
  
  // ---- no longer need the parcel; "Destroy" it. (clear() alone doesn't release the memory)
  Destroy(*this);  
  Destroy(kmer_locs_p);


}


























// ---------------------------------------------
// KmerParcelVec<K> (vec< KmerParcel<K> >)
// ---------------------------------------------

template<size_t K>
void KmerParcelVec<K>::Init(const size_t parcel0_ID,
                            const size_t parcel1_ID,
                            const size_t n_parcels,
                            const size_t n_blocks)
{
  //cout << "Init: " 
  //     << "parcel0_ID = " << parcel0_ID << " "
  //     << "parcel1_ID = " << parcel1_ID << " "
  //     << endl;

  _parcel0_ID = parcel0_ID;
  _parcel1_ID = parcel1_ID;
  _n_parcels = n_parcels;
  _n_blocks = n_blocks;

  _n_sub_parcels = parcel1_ID - parcel0_ID;
  (*this).resize(_n_sub_parcels);
    
  _n_sub_parcels_started = 0;
  _n_sub_parcels_done = 0;
    
  _n_blocks_started = 0;
  _n_blocks_done = 0;

}


template<size_t K>
void KmerParcelVec<K>::ParseReadKmersForParcelIDs(const BaseVecVec & bases,
                                                  const vec<size_t> & read1_IDs,
                                                  const size_t block_ID) 
{
  const size_t n_reads = bases.size();
  
  BaseVec kmer;
    

  // this block only looks at a read subset
  const size_t read0_ID = (block_ID == 0) ? 0 : read1_IDs[block_ID - 1];
  const size_t read1_ID = read1_IDs[block_ID];
  
  vec< vec< KmerRecord<K> > > parcel_bufs(_n_sub_parcels);


  // ---- estimate number of kmers for memory reservation; NOT TESTED so COMMENTED OUT FOR NOW
  //size_t n_kmer_estimate = ( 0.5 * (bases[read1_ID-1].size() + bases[read0_ID].size() - 2 * K) *
  //                           (read1_ID - read0_ID) ) / _n_parcels;

  //for (size_t sub_parcel_ID = 0; sub_parcel_ID != _n_sub_parcels; ++sub_parcel_ID)
  //  parcel_bufs[sub_parcel_ID].reserve(n_kmer_estimate);


  // ---- cycle through read subset
  for (size_t read_ID = read0_ID; read_ID != read1_ID; ++read_ID) {

    const BaseVec & current_read = bases[read_ID];
    if (current_read.size() >= K) {

      const size_t n_kmers = current_read.size() - K + 1;
        
      for (size_t base_ID = 0; base_ID < n_kmers; base_ID++) {
          
        // ribeiro: optimization opportunity: iain suggested not calculating 
        // the parcel_ID all the time.
        // It is enough to calculate a few bits of the parcel_ID to know if
        // the kmer belongs in the parcel or not; if it does, then we can
        // compute the remain bits.  
        const size_t parcel_ID = ComputeParcelIDFromBaseVec(current_read, base_ID, K, _n_parcels);
          
        if (parcel_ID >= _parcel0_ID && parcel_ID < _parcel1_ID) {
            
          size_t sub_parcel_ID = parcel_ID - _parcel0_ID;
            
          // All kmers are canonicalized currently, just to be sure.
          // ribeiro: There is a small optimization opportunity here, since, when 
          // computing the parcel_IDs, we get a clue of which kmers are canonical,
          // which ones are not, and which ones are undetermined.
          // However, it's not clear that this would improve performance significantly.
            
          const bool RC = (kmer.SetToSubOf(current_read, base_ID, K).Canonicalize() == BaseVec::rc_form);
            
          vec< KmerRecord<K> > & parcel_buf = parcel_bufs[sub_parcel_ID];
            
          // add kmer record to buffer
          parcel_buf.push_back(KmerRecord<K>(kmer, read_ID, base_ID, RC));
        
        }
      }
    }
  }
  // flush all buffers
  for (size_t sub_parcel_ID = 0; sub_parcel_ID != _n_sub_parcels; ++sub_parcel_ID) {
    (*this)[sub_parcel_ID].PushKmerRecords(parcel_bufs[sub_parcel_ID]);
    Destroy(parcel_bufs[sub_parcel_ID]); // really get rid of memory
  }
}




template<size_t K>
bool KmerParcelVec<K>::WaitingForBlocks()
{
  Locker lock(*this);
  return (_n_blocks_done < _n_blocks && _n_blocks_started == _n_blocks);
}


    

template<size_t K>
bool KmerParcelVec<K>::RunNextTask(const size_t thread_ID,
                                   const BaseVecVec & bases, 
                                   const vec<size_t> & read1_IDs,
                                   KmerParcelsStore & parcels_store,
                                   KmerFrequencyStatistics & stats)
{
  bool run_block = false;
  size_t block_ID = 0;
  bool run_parcel = false;
  size_t parcel_ID = 0;
  size_t sub_parcel_ID = 0;
  
  // 1st run blocks; all blocks must complete before starting sub_parcels
  // 2nd run sub parcels
  
  if (true) {
    Locker lock(*this);
    
    if (_n_blocks_done != _n_blocks) {

      if (_n_blocks_started != _n_blocks) {

        run_block = true;               // start a new block
        block_ID = _n_blocks_started;
        _n_blocks_started++;
      }

    }
    else { 

      if (_n_sub_parcels_started != _n_sub_parcels) {

        run_parcel = true;           // start a new parcel
        sub_parcel_ID = _n_sub_parcels_started;
        parcel_ID = _parcel0_ID + sub_parcel_ID;
        _n_sub_parcels_started++;
      }
    }

  } // data gets unlocked here.
 



  if (run_block) {
    //cout << Tag("RunNextTask") 
    //     << "running block " << block_ID << " "
    //     << "parcel01_IDs = " << _parcel0_ID << " " << _parcel1_ID << " "
    //     << endl;
    // parse a block of reads for kmers 
    ParseReadKmersForParcelIDs(bases, read1_IDs, block_ID);
    //cout << Tag("RunNextTask") 
    //     << "done with block " << block_ID << " "
    //     << "parcel01_IDs = " << _parcel0_ID << " " << _parcel1_ID << " "
    //     << endl;
      
    Locker lock(*this);
    _n_blocks_done++;
     
    return true;
  } 
    



  if (run_parcel) {
    //cout << Tag("RunNextTask") 
    //     << "sorting parcel " << parcel_ID << " "
    //     << "parcel01_IDs = " << _parcel0_ID << " " << _parcel1_ID << " " 
    //     << endl;
    (*this)[sub_parcel_ID].SortKmerParcel();
    (*this)[sub_parcel_ID].WriteSortedKmerParcel(thread_ID,
                                                 parcel_ID, 
                                                 parcels_store,
                                                 stats);
    //cout << Tag("RunNextTask") 
    //     << "done with parcel " << parcel_ID << " "
    //     << "parcel01_IDs = " << _parcel0_ID << " " << _parcel1_ID << " " 
    //    << endl;
    Locker lock(*this);
    _n_sub_parcels_done++;

    return true;
  }


  // nothing was done
  // either all parcels have already started 
  // or all blocks have started but haven't finished yet
  return false;
}
















// ---------------------------------------------
// KmerParcelVecVec<K> (vec< KmerParcelVec<K> >)
// ---------------------------------------------


template<size_t K>
KmerParcelVecVec<K>::KmerParcelVecVec(const size_t n_parcels,
                                      const size_t n_threads)
{
  size_t n_parcel_sets = (n_parcels + n_threads - 1) / n_threads;
  
  //cout << "n_sets = " << n_parcel_sets <<endl;
  //cout << "n_parcels = " << n_parcels << endl;
  //cout << "n_threads = " << n_threads << endl;
  
  (*this).resize(n_parcel_sets);
  for (size_t set_ID = 0; set_ID != n_parcel_sets; ++set_ID) {
    
    size_t parcel0_ID = set_ID * n_threads;
    size_t parcel1_ID = Min(parcel0_ID + n_threads, n_parcels);
    
    (*this)[set_ID].Init(parcel0_ID, 
                         parcel1_ID, 
                         n_parcels, 
                         n_threads);
  }

}






template<size_t K>
void KmerParcelVecVec<K>::RunTasks(const size_t thread_ID,
                                   const BaseVecVec & bases, 
                                   const vec<size_t> & read1_IDs,
                                   KmerParcelsStore & parcels_store,
                                   KmerFrequencyStatistics & stats)
                
{
  bool done = false;
    
  while (!done) {

    bool ran_something = false;
    typename vec< KmerParcelVec<K> >::iterator it = this->begin();
    bool ran;
    size_t set_ID = 0;
    do {
      ran = it->RunNextTask(thread_ID,
                            bases, 
                            read1_IDs,
                            parcels_store,
                            stats);
      //cout << Tag("RunTasks")
      //     << "thread_ID = " << thread_ID << " "
      //     << "set_ID = " << set_ID << " "
      //     << "ran = " << ran << endl;
        
      set_ID ++;
      it++;
        
    } while (!ran && it != this->end());
      
    if (!ran) {
      if (this->rbegin()->WaitingForBlocks()) sleep(5);
      else done = true;
    }

  }
}


















// ---------------------------------------------
// KmerParcels    Builder
// ---------------------------------------------




template<size_t K>
void KmerParcelsBuilder::BuildTemplate()
{
  //cout <<  "parcel_vec_vec<K>(" << _n_parcels << ", " << _n_threads << ")" << endl;

  KmerParcelVecVec<K> parcel_vec_vec(_n_parcels, _n_threads);
  
  ParcelProcessor<K> parcel_proc(_bases,
                                 _read1_IDs,
                                 &parcel_vec_vec,
                                 _parcels_store,
                                 _stats);
  
  if ( _n_threads <= 1 )
  {
      parcel_proc(0);
  }
  else
  {
      WorklistN< ParcelProcessor<K> > wl(parcel_proc,_n_threads,_n_threads-1);
  }
}












size_t KmerParcelsBuilder::EstimateTotalMemory()
{

  size_t n_bytes_kmer_record;
# define __SIZE(__K) n_bytes_kmer_record = sizeof(KmerRecord<__K>); 
  DISPATCH_ON_K( _K, __SIZE );

  size_t n_bytes = _n_kmers_total * n_bytes_kmer_record; 
  if (_verbose) {
    cout << Tag("EKPMB") << "n_reads = " << _bases.size() << endl;
    cout << Tag("EKPMB") << "n_kmers = " << _n_kmers_total << endl;
    cout << Tag("EKPMB") << "total memory = " 
         << static_cast<float>(n_bytes) / 1.0e9 << " GB" << endl;
  }
  return n_bytes;
}



// ---- Computes the number of kmer parcels so that the parcelization
//      only takes up to n_bytes_max extra memory.
void KmerParcelsBuilder::ComputeNumParcels()
{
  const size_t n_bytes = EstimateTotalMemory();
 
  // ---- minimum number of passes to process all kmers in parallel
  // the 2x factor below is necessary because, in the new (2010-03)
  // parallel scheme, data for up to 2 passes are stored in memory
  const size_t n_passes = (2 * n_bytes + _n_bytes_max - 1) / _n_bytes_max; 
  
  
  _n_parcels = _n_threads * n_passes;
  

  if (_verbose) {
    cout << Tag("KPB:CNKP") << "ComputeNumKmerParcels()" << endl;
    cout << Tag("KPB:CNKP") << "n_reads = " << _bases.size() << endl;
    cout << Tag("KPB:CNKP") << "K = " << _K << endl;
    cout << Tag("KPB:CNKP") << "n_GB_max = " << static_cast<double>(_n_bytes_max) / (1 << 30) << endl;
    cout << Tag("KPB:CNKP") << "n_GB_kmer_parcels = " << static_cast<double>(n_bytes) / (1 << 30) << endl;
    cout << Tag("KPB:CNKP") << "n_parcels = " << _n_parcels << endl;
    cout << Tag("KPB:CNKP") << "n_passes = " << n_passes << endl;
  }

  if (_K < 16) {
    // make sure K is large enough to generate parcel_IDs
    // K = 2,3 => np_max =   8 = 1 << 3
    // K = 4,5 => np_max =  64 = 1 << 6
    // K = 6,7 => np_max = 512 = 1 << 9
    size_t n_parcels_max = 1ul << (3 * (_K >> 1));
    ForceAssertLe(_n_parcels, n_parcels_max);  
  }
}
 






void  KmerParcelsBuilder::ComputeRead1IDs()
{

  const size_t n_reads = _bases.size();
  _read1_IDs.resize(_n_threads);

  /*
  cout << "n_reads = " << n_reads << endl;
  cout << "n_kmers_total = " << _n_kmers_total << endl;
  cout << "n_threads = " << _n_threads << endl;
  */

  size_t n_kmers = 0;
  size_t read_ID = 0;
  for (size_t block_ID = 0; block_ID < _n_threads; block_ID++) {

    const size_t n_kmers_cum = (_n_kmers_total * (block_ID + 1)) / _n_threads;

    while (read_ID != n_reads && n_kmers < n_kmers_cum) {
      const size_t n_bases = _bases[read_ID++].size();
      if (n_bases >= _K) 
        n_kmers += n_bases - _K + 1;
    }
    _read1_IDs[block_ID] = read_ID;
  }

  /*
  cout << "n_blocks = " << _read1_IDs.size() << endl;
  for (size_t block_ID = 0; block_ID < _n_threads; block_ID++) {
    cout << "read1_IDs[" << block_ID << "] = " << _read1_IDs[block_ID] 
         << " (" << (_read1_IDs[block_ID] - ((block_ID == 0) ? 0 : _read1_IDs[block_ID - 1])) << ")"
         << endl;
  }
  */
}






void KmerParcelsBuilder::Build(const size_t n_parcels)
{
  const size_t n_reads = _bases.size();
  if (n_reads == 0) {
    cout << Tag("KPB:B") << "NOTHING TO DO: n_reads = 0." << endl;
    ForceAssertGt(n_reads, 0u);
  }

  _n_parcels = n_parcels;

  if (_n_parcels == 0) 
    ComputeNumParcels();

  if (_verbose) {
    cout << Tag("KPB:B") << "start" << endl;
    cout << Tag("KPB:B") << "GREPABLE: K = " << _K << endl;
    cout << Tag("KPB:B") << "GREPABLE: num_reads = " << n_reads << endl;
    cout << Tag("KPB:B") << "GREPABLE: num_kmers = " << _n_kmers_total << endl;
    cout << Tag("KPB:B") << "GREPABLE: num_parcels = " << _n_parcels << endl;
    cout << Tag("KPB:B") << "n_threads = " << _n_threads << endl;
  }
  
  if (_n_parcels == 0) {
    if (_verbose)
      cout << Tag("KPB:B") << "n_parcels = 0, nothing to do." << endl;
  }
  else {
  
    if (_parcels_store.GetNumParcels() != 0) 
      cout << Tag("KPB:B") << "WARNING: parcels already exist! We are about to remake them!!!" << endl;
    
    
    ComputeRead1IDs();
  

    // ribeiro: Templatized in K solely because of KmerRecord<K>.
    // ribeiro: sorting a vec< KmerRecord<K> > is 2 to 3 times faster than
    // ribeiro: sorting a vec< pair<BaseVec, KmerLoc> >

  
    _parcels_store.SetNumParcels(_n_parcels);
    _parcels_store.Open();
  
# define __BKT(__K) BuildTemplate<__K>()

    DISPATCH_ON_K( _K, __BKT );
  
    _parcels_store.Close();
  





    // ---- collecting kmer frequency stats
  
    String dir_name = _parcels_store.GetDirectoryName();
  
    if (_verbose)
      cout << Tag("KPB:B") << "collecting reads min and max kmer frequency data" << endl;
  
    if (_verbose || dir_name != "") {
      MapOfCounters reads_min_kmer_freq_count(_stats.reads_min_kmer_freq);
      MapOfCounters reads_max_kmer_freq_count(_stats.reads_max_kmer_freq);
    
      /*
      if (dir_name != "") {
        reads_min_kmer_freq_count.Write(dir_name + "/reads_min_kmer_freq.count.dat");
        reads_max_kmer_freq_count.Write(dir_name + "/reads_max_kmer_freq.count.dat");
      }
      */
    
      if (_verbose) {
        PrintValueAndPercentage(cout, Tag("KPB:B") + "GREPABLE: num_reads_min_kmer_freq_1 = ",
                                reads_min_kmer_freq_count[1], n_reads);
        PrintValueAndPercentage(cout, Tag("KPB:B") + "GREPABLE: num_reads_min_kmer_freq_2 = ",
                                reads_min_kmer_freq_count[2], n_reads);
        PrintValueAndPercentage(cout, Tag("KPB:B") + "GREPABLE: num_reads_max_kmer_freq_1 = ",
                                reads_max_kmer_freq_count[1], n_reads);
        PrintValueAndPercentage(cout, Tag("KPB:B") + "GREPABLE: num_reads_max_kmer_freq_2 = ",
                                reads_max_kmer_freq_count[2], n_reads);
      }
    }
  
  
    if (_verbose) {
      cout << Tag("KPB:B") << "collecting kmer frequency data" << endl;
      size_t n_distinct_kmers = _stats.AllKmerFrequencyCounters().Total();
      PrintValueAndPercentage(cout, Tag("KPB:B") + "GREPABLE: num_distinct_kmers        = ",
                              n_distinct_kmers, n_distinct_kmers);
      PrintValueAndPercentage(cout, Tag("KPB:B") + "GREPABLE: num_distinct_kmers_freq_1 = ",
                              _stats.AllKmerFrequencyCounters()[1], n_distinct_kmers);
      PrintValueAndPercentage(cout, Tag("KPB:B") + "GREPABLE: num_distinct_kmers_freq_2 = ",
                              _stats.AllKmerFrequencyCounters()[2], n_distinct_kmers);
    }
  
    /*
    if (dir_name != "") {
      _stats.AllKmerFrequencyCounters().Write(dir_name + "/kmer_frequencies.count.dat");
      const vec<KmerFrequencyCounters> & gc_kfc = _stats.AllKmerGCFrequencyCounters();
      
      if (_write_gc_stats) {
        const size_t n = gc_kfc.size();
        for (size_t i = 0; i != n; i++) {
          gc_kfc[i].Write(dir_name + "/kmer_frequencies.count.gc=" + ToString(i) + ".dat");
        } 
      }
    }
    */
  
    if (_verbose)
      cout << Tag("KPB:B") << "end" << endl;

  }
}




