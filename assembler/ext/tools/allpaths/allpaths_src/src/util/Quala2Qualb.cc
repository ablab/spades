// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research

// Quala2Qualbb: convert a fasta quality score file into a qualb file.
// Like Qualb, but without the need for PRE or a special extension.

#include <ctype.h>

#include <strstream>

#include "MainTools.h"
#include "FastIfstream.h"
#include "Qualvector.h"

int main(int argc, char *argv[])
{
  RunTime();

  BeginCommandArguments;
  CommandArgument_String(IN);
  CommandArgument_String(OUT);

  // Read trimming and reversing
  CommandArgument_Int_OrDefault(SKIP_RIGHT_BASES, 0);
  CommandArgument_Int_OrDefault(SKIP_LEFT_BASES, 0);
  CommandArgument_Int_OrDefault_Doc(END_BASE, 0, "If 0, don't trim");
  CommandArgument_Bool_OrDefault( REVERSE_READS, False);
  EndCommandArguments;


  QualVecVec quals_out;

  int total_seqs = 0;
  longlong total_bases = 0;
  for (int pass = 1; pass <= 2; pass++) {
    if (pass == 2) quals_out.Reserve(total_bases, total_seqs);
    fast_ifstream quals_in(IN);
    String line;
    QualVec q;
    Bool first = True;
    while(1) {
      getline(quals_in, line);
      if (quals_in.fail()) {
        if (pass == 1) {
          total_bases += q.size();
          ++total_seqs;
        }
        if (pass == 2) 
          quals_out.push_back(q);
        q.clear();
        break;
      }
      ForceAssert(line.size() > 0);
      if (line[0] == '>') {
        if (!first) {
          if (pass == 1) {
            total_bases += q.size();
            ++total_seqs;
          }
          if (pass == 2) 
            quals_out.push_back(q);
        }
        first = False;
        q.clear();
        while(1) {
          char c;
          quals_in.peek(c);
          if (quals_in.fail() || c == '>') break;
          getline(quals_in, line);
          for (int j = 0; j < (int) line.size(); j++)
            ForceAssert(isspace(line[j]) || isdigit(line[j]));
          istrstream i(line.c_str());
          while(1) {
            int n;
            i >> n;
            if (i.fail()) break;
            q.push_back(n);
          }
        }
      }
      if (quals_in.fail()) {
        if (pass == 1) {
          total_bases += q.size();
          ++total_seqs;
        }
        if (pass == 2) 
          quals_out.push_back(q);
        q.clear();
        break;
      }
    }
  }


  if ((SKIP_RIGHT_BASES != 0) || 
      (SKIP_LEFT_BASES != 0) ||
      (END_BASE != 0))
  {
    for (size_t i = 0; i < quals_out.size(); i++) {

        int start = SKIP_LEFT_BASES;
        int len   = quals_out[i].size() - SKIP_RIGHT_BASES - SKIP_LEFT_BASES;
        if (END_BASE) {
          len = END_BASE + 1 - start;
          if (len + start > (signed)quals_out[i].size())
            len = quals_out[i].size() - start;
        }
            
        quals_out[i].SetToSubOf(quals_out[i], start, len);

        if (REVERSE_READS) 
          quals_out[i].ReverseMe();

      }
  }


          
  quals_out.WriteAll(OUT);
}
