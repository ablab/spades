///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/**
 * FastbQualbTrimReverse
 *
 * IN_HEAD       it needs <IN_HEAD>.{fastb,qualb}
 * OUT_HEAD      it saves <OUT_HEAD>.{fastb,qualb}
 * TRIM_START    first base to keep
 * TRIM_END      last base to keep
 * REVERSE       wether to reverse the reads or not
 *
 * \author Filipe Ribeiro    11-2010
 */

#include "MainTools.h"
#include "Basevector.h"
#include "Qualvector.h"


static inline 
String Tag(String S = "FQTR") { return Date() + " (" + S + "): "; } 

int main(int argc, char *argv[])
{
  RunTime();
  
  BeginCommandArguments;
  CommandArgument_String(IN_HEAD);
  CommandArgument_String(OUT_HEAD);
  CommandArgument_Int_OrDefault(TRIM_START, 0);
  CommandArgument_Int_OrDefault(TRIM_END, 0);
  CommandArgument_Bool_OrDefault(REVERSE, False);
  EndCommandArguments;
  
  // Dir and file names.
  String fn_reads_in = IN_HEAD + ".fastb";
  if (! IsRegularFile(fn_reads_in)) {
    cout << Tag() << "'" << fn_reads_in << "' not found. Abort." << endl;
    return 0;
  }

  bool do_quals = true;
  String fn_quals_in = IN_HEAD + ".qualb";
  if (! IsRegularFile(fn_quals_in)) {
    cout << Tag() << "'" << fn_quals_in << "' not found." << endl;
    do_quals = false;
  }
  String out_dir = ".";
  if (OUT_HEAD.Contains("/")) out_dir = OUT_HEAD.RevBefore("/");
  Mkpath(out_dir);
  
  String fn_bases_out = OUT_HEAD + ".fastb";
  String fn_quals_out = OUT_HEAD + ".qualb";
  

  size_t n_reads;
  vec<size_t> lens_reads;

  // ---- bases
  {
    cout << Tag() << "Trimming bases to '" << fn_bases_out << "'." << endl;
    BaseVecVec bases(fn_reads_in);
    IncrementalWriter<BaseVec> bwriter(fn_bases_out.c_str());
    n_reads = bases.size();
    for (size_t i = 0; i < n_reads; i++) {
      if (do_quals)
        lens_reads.push_back(bases[i].size());

      const int end = ((TRIM_END > 0 && TRIM_END < int(bases[i].size())) ? 
                       TRIM_END : bases[i].size() - 1);
      const int len = end - TRIM_START + 1;

      bases[i].SetToSubOf(bases[i], TRIM_START, len);

      if (REVERSE)
        bases[i].ReverseComplement();
      
      bwriter.add(bases[i]);
    }
    bwriter.close();
  }

  // ---- quals
  if (do_quals) {
    size_t n_warns = 0;
    cout << Tag() << "Trimming quals to '" << fn_quals_out << "'." << endl;
    QualVecVec quals(fn_quals_in);
    IncrementalWriter<QualVec> qwriter(fn_quals_out.c_str());
    for (size_t i = 0; i < n_reads; i++) {
      
      if (lens_reads[i] != quals[i].size()) {
        if (n_warns < 20) {
          cout << Tag() << "WARNING: inconsistent sizes for read " << i << ": " 
               << lens_reads[i] << " bases and " << quals[i].size() << " quals." << endl; 
        }
        n_warns++;        
      }

      const int end = ((TRIM_END > 0 && TRIM_END < int(quals[i].size())) ? 
                       TRIM_END : quals[i].size() - 1);

      const int len = end - TRIM_START + 1;

      quals[i].SetToSubOf(quals[i], TRIM_START, len);

      if (REVERSE)
        quals[i].ReverseMe();
    
      qwriter.add(quals[i]);
    }
    qwriter.close();

    if (n_warns > 0)
      cout << Tag() << "WARNING: found " << n_warns << " inconsistent read/qual sizes." << endl;

  }  
  

  // ---- done
  cout << Tag() << "Done." << endl;
  
}

