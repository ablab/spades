/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////



/// Transform one fasta file into fastb format, save names separately.
/// 
/// \file Fasta2Fastb.cc
///
/// names are saved as a vecString in OUT.names



#include "Basevector.h"
#include "MainTools.h"
#include "FastaFileset.h"

int main(int argc, char *argv[])
{
  RunTime();

  BeginCommandArguments;
  CommandArgument_String(IN);
  CommandArgument_String_OrDefault(OUT,"");

  // Read trimming and reversing
  CommandArgument_Int_OrDefault(SKIP_RIGHT_BASES, 0);
  CommandArgument_Int_OrDefault(SKIP_LEFT_BASES, 0);
  CommandArgument_Int_OrDefault_Doc(END_BASE, 0, "If 0, don't trim");
  CommandArgument_Bool_OrDefault( REVERSE_READS, False);

  CommandArgument_Bool_OrDefault_Doc(NAMES, True,
		     "Whether to save the original sequence names as found in "
				     "fasta file into fastb.names file");
  EndCommandArguments;

  BaseVecVec reads;
  vecString readNames;

  if (OUT.empty()) OUT=IN.SafeBefore(".fa") + ".fastb";

  FastFetchReads(reads, &readNames, IN);
  if ((SKIP_RIGHT_BASES != 0) || 
      (SKIP_LEFT_BASES != 0) ||
      (END_BASE != 0))
  {
    for (size_t i = 0; i < reads.size(); i++) {

        int start = SKIP_LEFT_BASES;
        int len   = reads[i].size() - SKIP_RIGHT_BASES - SKIP_LEFT_BASES;
        if (END_BASE) {
          len = END_BASE + 1 - start;
          if (len + start > (signed)reads[i].size())
            len = reads[i].size() - start;
        }
            
        reads[i].SetToSubOf(reads[i], start, len);

        if (REVERSE_READS) 
          reads[i].ReverseComplement();

      }
  }
  reads.WriteAll(OUT);
  if (NAMES) readNames.WriteAll(OUT + ".names");
  return 0;
}
