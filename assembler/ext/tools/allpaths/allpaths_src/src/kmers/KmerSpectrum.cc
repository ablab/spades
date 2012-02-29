///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#include "MainTools.h"

#include "Basevector.h"
#include "Qualvector.h"
#include "feudal/VirtualMasterVec.h"

static inline 
String Tag(String S = "KS") { return Date() + " (" + S + "): "; } 

#include "kmers/KmerSpectra.h"

#include "kmers/KmerSpectrumCore.h"


typedef VirtualMasterVec<BaseVec> VBaseVecVec;


size_t number_of_bases_sample(const String fn)
{
  VBaseVecVec vbvv(fn.c_str()); 
  const size_t nr = vbvv.size();
  const size_t ns = 500;
  if (ns > nr) {
    size_t nb = 0;
    for (VBaseVecVec::const_iterator it = vbvv.begin(); 
         it != vbvv.end(); it++)
      nb += it->size();
    return nb;
  }
  else {
    const size_t dir = nr / ns;
    const size_t ir0 = dir / 2;
    size_t nbs = 0;
    for (size_t ir = ir0; ir < nr; ir += dir)
      nbs += vbvv[ir].size();

    return (nr * nbs) / ns;
  }
  

}







int main(int argc, char *argv[])
{
  RunTime();

  BeginCommandArguments;
  CommandArgument_UnsignedInt_Doc(K, "Kmer size.");
  CommandArgument_String_Doc(HEAD, "looks for '<HEAD>.fastb'.");
  CommandArgument_String_OrDefault_Doc(HEAD2, "", "looks for '<HEAD2>.fastb'.");
  CommandArgument_Bool_OrDefault_Doc(DO_QUALS, False, "Whether to separate the spectrum by kmer quality."); 
  CommandArgument_Bool_OrDefault_Doc(DO_AFFIXES, False, "Whether to separate the spectrum by kmer affixes."); 

  CommandArgument_UnsignedInt_OrDefault_Doc(KF_MIN, 1, "Minimum kmer frequency when computing affixes spectra.");
  CommandArgument_UnsignedInt_OrDefault_Doc(KF_MAX, 0, "Maximum kmer frequency when computing affixes spectra.");
  
  CommandArgument_Bool_OrDefault_Doc(G_ESTIMATE, False, "When true compute genome size estimate.");
  CommandArgument_Bool_OrDefault_Doc(G_ESTIMATE_ONLY, False, "When true don't compute spectrum but load it from '<HEAD>.<K>mer.kspec' and compute genome size estimate.");
  CommandArgument_UnsignedInt_OrDefault_Doc(KF_LOW, 10, "Under-estimate of kmer frequency for local kmer spectrum minimum (for genome size estimate).");
  CommandArgument_UnsignedInt_OrDefault_Doc(READ_LEN, 0, "Optional read length in data (for genome size estimate).");

  CommandArgument_UnsignedInt_OrDefault_Doc(PLOIDY, 1, "If 2 will estimate SNP rate.");
  CommandArgument_String_OrDefault_Doc(PLOIDY_FILE, "", "If file exists reads ploidy and, if ploidy = 2, estimates SNP rate.");

  CommandArgument_UnsignedInt_OrDefault_Doc(VERBOSITY, 1, "Level of verbosity.");

  CommandArgument_LongLong_OrDefault(MAX_MEMORY_GB, 0);
  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, "Number of threads to use in parallel");
  EndCommandArguments;

  // ---- Thread control

  SetMaxMemory(MAX_MEMORY_GB << 30);
  NUM_THREADS = configNumThreads(NUM_THREADS);

  // ---- Ploidy

  if (PLOIDY_FILE != "" && IsRegularFile(PLOIDY_FILE)) // superseeds PLOIDY argument
    PLOIDY = StringOfFile(PLOIDY_FILE, 1).Int();

  const String fn_fastb = HEAD + ".fastb";
  
  // ----  do kmer spectrum on one file only

  if (HEAD2 == "") {
  
      // ---- do kmer spectra by qualities

    if (DO_QUALS) { 

      cout << Tag() << "Loading bases from '" << fn_fastb << "'." << endl;
      BaseVecVec bvv(fn_fastb);
      
      const String fn_qualb = HEAD + ".qualb";
      cout << Tag() << "Loading quals from ." << fn_qualb << "'." << endl;
      QualVecVec qvv(fn_qualb);
      
      KmerQualitySpectra kqspec(K);
      kmer_spectra_by_quals_compute(bvv, qvv, &kqspec, VERBOSITY, NUM_THREADS);
      
      cout << Tag() << "Writing kmer quality spectra." << endl;
      kqspec.to_text_file(HEAD);
    }

    // ---- don't do quals, but do affixes
    
    else if (DO_AFFIXES) {

      if (KF_MAX > 1) ForceAssertLe(KF_MIN, KF_MAX);

      KmerSpectraAffixes kspecs(K);
      Validator validator(KF_MIN, KF_MAX);

      // ---- read data

      cout << Tag() << "Loading bases from '" << fn_fastb << "'." << endl;
      BaseVecVec bvv(HEAD + ".fastb");

      // ---- compute affixes spectra

      cout << Tag() << "Computing kmer spectra by affixes for kmer frequency range"
           << " [" << KF_MIN << ":";
      if (KF_MAX == 0) cout << "inf.]." << endl;
      else             cout << KF_MAX << "]." << endl;
      
      kmer_spectra_by_affixes_compute(bvv, & kspecs, validator, VERBOSITY, NUM_THREADS);


      // ---- write kmer spectra
      
      cout << Tag() << "Writing kmer affixes spectra." << endl;
      
      {
        const KmerSpectrum kspec_all = kspecs.all();
        kspec_all.to_text_file(HEAD);
        cout << Tag() << "nk(all)= " << kspec_all.sum() << endl;
      }
      for (unsigned n_pre = 0; n_pre <= 4; n_pre++) {
        for (unsigned n_suf = n_pre; n_suf <= 4; n_suf++) {
          const KmerSpectrum & kspec = kspecs(n_pre, n_suf);
          kspec.to_text_file(HEAD + "." + ToString(n_pre) + "-" + ToString(n_suf));
          cout << Tag() << "nk(" << n_pre << "-" << n_suf << ")= " << kspec.sum() << endl;
        }
      }

    }

    // ---- do just plain kmer spectrum

    else {
      KmerSpectrum kspec(K);
      size_t read_len = 0;

      if (G_ESTIMATE_ONLY) {
        
        // ---- get number of bases

        if (READ_LEN) {
          read_len = READ_LEN;
        } 
        else {
          cout << Tag() << "Counting reads in '" << fn_fastb << "'." << endl;
          size_t nr = MastervecFileObjectCount(fn_fastb);
          cout << Tag() << "Estimating number of bases." << endl;
          size_t nb = number_of_bases_sample(fn_fastb);
          read_len = nb / nr;
        }

        // ---- read kmer spectrum from text file
        cout << Tag() << "Reading kmer spectrum." << endl;
        kspec.from_text_file(HEAD);

      }
      else {
        // ---- read data
        cout << Tag() << "Loading bases from '" << fn_fastb << "'." << endl;
        BaseVecVec bvv(fn_fastb);
        
        // ---- compute spectrum
        cout << Tag() << "Computing kmer spectrum." << endl;
        kmer_spectrum_compute(bvv, & kspec, VERBOSITY, NUM_THREADS);
        
        cout << Tag() << "Writing kmer spectrum." << endl;
        kspec.to_text_file(HEAD);
        
        // ---- compute number of bases
        size_t nr = bvv.size();
        size_t nb = 0;
        for (size_t ir = 0; ir != nr; ir++) 
          nb += bvv[ir].size();
        read_len = nb / nr;
      }

      // ---- estimate genome size
      
      if (G_ESTIMATE || G_ESTIMATE_ONLY)
        genome_analysis_report(kspec, read_len, PLOIDY, KF_LOW, VERBOSITY);
      
    }
  }


  // ----  do kmer spectra of two data sets together 

  else {
    cout << Tag() << "Loading first set of bases from '" << fn_fastb << "'." << endl;
    BaseVecVec bvv(fn_fastb);
    const size_t nbv1 = bvv.size();

    const String fn2_fastb = HEAD2 + ".fastb";
    cout << Tag() << "Loading second set of bases from '" << fn2_fastb << "'." << endl;
    VBaseVecVec vbvv2(fn2_fastb.c_str());
    bvv.append(vbvv2.begin(), vbvv2.end());

  
    cout << Tag() << "Computing kmer bi-spectrum." << endl;
        
    KmerBiSpectrum kbspec(K);
    kmer_bi_spectrum_compute(bvv, nbv1, &kbspec, VERBOSITY, NUM_THREADS);
        
    cout << Tag() << "Writing kmer bispectra." << endl;
  
    kbspec.AB.to_text_file(HEAD + "." + HEAD2);

    kbspec.A_in.to_text_file(HEAD + ".in." + HEAD2);
    kbspec.A_out.to_text_file(HEAD + ".out." + HEAD2);
  
    KmerSpectrum kspecA = kbspec.A_in + kbspec.A_out;
    kspecA.to_text_file(HEAD);
    
    kbspec.B_in.to_text_file(HEAD2 + ".in." + HEAD);
    kbspec.B_out.to_text_file(HEAD2 + ".out." + HEAD);
    
    KmerSpectrum kspecB = kbspec.B_in + kbspec.B_out;
    kspecB.to_text_file(HEAD2);
  }
  
}
