/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/// Transform fastb file IN into fasta format.
/// 
/// \file Fastb2Fasta.cc
///
/// Default for OUT is IN's prefix + .fasta
/// If GENERATE_NAMES is set names must be in file IN.names or IN.names.gz.

#include "Basevector.h"
#include "MainTools.h"

int main(int argc, char *argv[])
{
  RunTime();

  BeginCommandArguments;
  CommandArgument_String(IN);
  CommandArgument_String_OrDefault_Doc(OUT, "", "use '-' for stdout.");
  CommandArgument_Bool_OrDefault(GENERATE_NAMES, True);
  CommandArgument_Bool_OrDefault(ONE_LINE_PER_SEQUENCE, False);
  CommandArgument_Int_OrDefault(BREAKCOL, 80);
  CommandArgument_UnsignedInt_OrDefault_Doc(READ_COUNT, 0,
    "Only convert READ_COUNT basevectors then stop.");
  EndCommandArguments;

  if (OUT.empty()) 
    OUT = IN.SafeBefore(".fastb") + ".fasta";

  // Setup virtual master vec to stream in the fastb
  typedef VirtualMasterVec<basevector> VmvBv_t;
  VmvBv_t reads( IN.c_str() );

  // Load 
  vecString readNames;
  if (!GENERATE_NAMES ) {
    readNames.ReadAll(IN + ".names");
  }

  const bool to_stdout = (OUT == "-");

  size_t count = 0;
  ofstream ofs;
  if (!to_stdout) ofs.open(OUT.c_str());
  ostream & os = (to_stdout) ? cout : ofs;

  for (VmvBv_t::Itr reads_itr = reads.begin(); reads_itr != reads.end() ; ++reads_itr) {
    if (ONE_LINE_PER_SEQUENCE) {
      os << reads_itr->ToString() << endl;
      count++;
    } 
    else {
      reads_itr->PrintCol(os, ToString(count++), BREAKCOL);
    }
    if (count == READ_COUNT)
      break;
  }

  if (!to_stdout) ofs.close();
}
