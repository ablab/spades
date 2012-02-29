///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// FastbQualbToFastq.  Given fastb and qualb, generate fastq file(s). 

#include "MainTools.h"
#include "Basevector.h"
#include "Qualvector.h"


String qual_vec_to_fastq_string(const QualVec & qv, const unsigned ascii_offset)
{
  ForceAssert(ascii_offset == 33 || ascii_offset == 64);
  String s = "";
  for (size_t i = 0; i < qv.size(); i++) 
    s += qv[i] + ascii_offset;
  return s;
}




void read_to_fastq_stream(ofstream & fastq,
                          const BaseVec & bv,
                          const QualVec & qv,
                          const size_t ibv,
                          const unsigned phred_offset,
                          const String & prefix,
                          const String & suffix,
                          bool rc)
{
  fastq << "@" << prefix << ":" << ibv << suffix << endl;
  fastq << (rc ? BaseVec(bv).ReverseComplement() : bv).ToString() << endl;
  fastq << "+" << endl;
  fastq << qual_vec_to_fastq_string((rc ? QualVec(qv).ReverseMe() : qv),
                                      phred_offset) << endl;
}                          






typedef VirtualMasterVec<BaseVec> VBaseVecVec;
typedef VirtualMasterVec<QualVec> VQualVecVec;

int main(int argc, char *argv[])
{
  RunTime();

  BeginCommandArguments;
  CommandArgument_String_Doc(HEAD_IN, "Need HEAD_IN.{fastb,qualb}");
  CommandArgument_String_Doc(HEAD_OUT, "Create HEAD_OUT.fastq");
  CommandArgument_Bool_Doc(PAIRED, "Whether the reads in the fastb/qualb are paired (interleaved).");
  CommandArgument_UnsignedInt_Doc(PHRED_OFFSET, "Either 33 or 64");
  CommandArgument_String_OrDefault_Doc(NAMING_PREFIX, HEAD_IN, "Read name prefix to read ID");
  CommandArgument_String_OrDefault_Doc(NAMING_SUFFIX, "", "Read name suffix to read ID");
  CommandArgument_Bool_OrDefault_Doc(PICARD_NAMING_SCHEME, False, "If true, "
				     "add '/1' to the end of first read "
				     "names, '/2' to the end of second.");
  CommandArgument_Bool_OrDefault_Doc(FLIP, False, "reverse complement reads and quals");
  EndCommandArguments;

  // Load the data.

  ForceAssert(PHRED_OFFSET == 33 || PHRED_OFFSET == 64);

  VBaseVecVec vbvv((HEAD_IN + ".fastb").c_str());
  VQualVecVec vqvv((HEAD_IN + ".qualb").c_str());

  VBaseVecVec::const_iterator it_bv = vbvv.begin();
  VQualVecVec::const_iterator it_qv = vqvv.begin();
    
  const size_t nbv = vbvv.size();
  size_t ibv = 0;

  if (PAIRED) {

    cout << "Creating '" << HEAD_OUT << ".{A,B}.fastq':" << endl;
    ofstream fastq_A((HEAD_OUT + ".A.fastq").c_str());
    ofstream fastq_B((HEAD_OUT + ".B.fastq").c_str());
    
    const String suffix_A = (PICARD_NAMING_SCHEME ? "/1" : "" ) + NAMING_SUFFIX;
    const String suffix_B = (PICARD_NAMING_SCHEME ? "/2" : "" ) + NAMING_SUFFIX;

    while (it_bv != vbvv.end() && it_qv != vqvv.end()) {
      
      ForceAssertEq(it_bv->size(), it_qv->size());
      
      read_to_fastq_stream(fastq_A, *it_bv, *it_qv, ibv, 
                           PHRED_OFFSET, NAMING_PREFIX, suffix_A, FLIP);
      it_bv++;
      it_qv++;
      
      read_to_fastq_stream(fastq_B, *it_bv, *it_qv, ibv, 
                           PHRED_OFFSET, NAMING_PREFIX, suffix_B, FLIP);
      it_bv++;
      it_qv++;

      dots_pct(ibv++, nbv/2);
    }

    fastq_A.close();
    fastq_B.close();

  }
  else {
    cout << "Creating '" << HEAD_OUT << ".fastq':" << endl;
    ofstream fastq((HEAD_OUT + ".fastq").c_str());
    
    while (it_bv != vbvv.end() && it_qv != vqvv.end()) {
      
      ForceAssertEq(it_bv->size(), it_qv->size());

      read_to_fastq_stream(fastq, *it_bv, *it_qv, ibv, 
                           PHRED_OFFSET, NAMING_PREFIX, NAMING_SUFFIX, FLIP);
      it_bv++;
      it_qv++;

      dots_pct(ibv++, nbv);
    }
    fastq.close();
  }

}
