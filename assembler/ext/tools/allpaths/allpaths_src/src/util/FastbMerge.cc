///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


/* FastbMerge
 *
 * Takes a set of fastb files and merges them into one. 
 *
 * Optionally merges quality scores in corresponding qualb files.
 *
 *
 * 2011-07    Filipe Ribeiro     <crdhelp@broadinstitute.org>
 *
 ******************************************************************************/


#include "MainTools.h"
#include "Basevector.h"
#include "Qualvector.h"
#include "String.h"

static inline 
String Tag(String S = "FM") { return Date() + " (" + S + "): "; } 







int main(int argc, char **argv)
{
  RunTime();
  
  BeginCommandArguments;
  CommandArgument_StringSet(HEADS_IN);
  CommandArgument_String(HEAD_OUT);
  EndCommandArguments;

  const size_t n_files = HEADS_IN.size();

  vec<String> fns_bases_in;
  vec<String> fns_quals_in;

  size_t n_reads_out = 0;
  bool do_quals = true;
  cout << Tag() << "Data to merge:" << endl;
  for (size_t i = 0; i < n_files; i++) {
    const String head_in = HEADS_IN[i];
    fns_bases_in.push_back(head_in + ".fastb");
    fns_quals_in.push_back(head_in + ".qualb");

    const size_t n_r = MastervecFileObjectCount(fns_bases_in[i]);
    cout << Tag() << setw(14) << n_r << "  reads in '" << fns_bases_in[i] << "'." << endl;
   
    if (IsRegularFile(fns_quals_in[i]))
      ForceAssertEq(n_r, MastervecFileObjectCount(fns_quals_in[i]));
    else 
      do_quals = false;

    n_reads_out += n_r;
  }
  cout << Tag() << setw(14) << n_reads_out << "  total reads." << endl;

  const String fn_bases_out = HEAD_OUT + ".fastb";
  const String fn_quals_out = HEAD_OUT + ".qualb";

  cout << Tag() << "Merging read set into '" << fn_bases_out << "'." << endl;
  MergeMastervecs(fn_bases_out, fns_bases_in);

  if (do_quals) {
    cout << Tag() << "Merging qual set into '" << fn_quals_out << "'." << endl;
    MergeMastervecs(fn_quals_out, fns_quals_in);
  }

  
  cout << Tag() << "Done." << endl;

}
