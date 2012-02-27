/***************************************************************************
 * Title:          SeqUtils.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/28/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "SeqUtils.h"
#include "utils.h"
#include "SimpleStats.h"

void PrintDNASequences(DNASequenceList &sequences, std::string &outFileName) {
	std::ofstream outFile;
	openck(outFileName, outFile, std::ios::out);
	ssize_t i;
	for (i = 0; i < sequences.size(); i++) {
		sequences[i].PrintSeq(outFile);
	}
}



char RandomNuc() {
	return nuc_char[Random(4)];
}


char MutateNuc(char c) {
  char m = c;
  while (m == c) {
    m = RandomNuc();
  }
	m = tolower(m);
  return m;
}


void MakeRC(DNASequence &src, DNASequence &rev) {
  unsigned char *revSeq;
  rev.strand = DNASequence::REVERSE_STRAND;
  if (src._ascii) {
    MakeRC((char*) src.seq, src.length, revSeq);
    rev._ascii = 1;
  }
  else {
    // Need to make a binary reverse complement.
    ssize_t s;
    char n, c;
    revSeq = new unsigned char[src.length];
    for (s = 0; s < src.length; s++) {
      //      n = DNASequence::UnmaskedValue(src.seq[src.length-s-1]);
      n = src.seq[src.length-s-1];
      if (n < -12 || n > 15) {
				std::cout << "error making reverse complement " << std::endl;
				exit(0);
      }
      c = comp_bin[(unsigned char) n];
      if (c != 16)
				revSeq[s] = c;
      else {
				std::cout  << "nonstandard nucleotide in rc, exiting " << std::endl;
				exit(0);
      }
    }
  }
  rev.seq = (unsigned char*) revSeq;
  rev.length = src.length;
}

  

char NucRC(char nuc) {
  return comp_bin[(unsigned char) numeric_nuc_index[(unsigned char) nuc]];
}

// Reverse complement code for ascii sequences
// TODO: Why "char *" for src and "unsigned char*" for dst?  The calls to this are really messy because of that.
void MakeRC(char* seq, ssize_t len, unsigned char* &rev){
  ssize_t s;
  unsigned char n, c;
  //UNUSED// char *complement;
  // Allocate a new sequence.
  rev = new unsigned char[len+1];
  for (s = 0; s < len; s++) {
    n = seq[len-s-1];
    c = comp_ascii[n];
    if (c != 0)
      rev[s] = c;
    else 
      // Nonstandard nuleotide (masked, etc). Leave as is.
      rev[s] = seq[len-s-1];
  }
}

void CountACTG(DNASequence &seq, ssize_t &a, ssize_t &c, ssize_t &t, ssize_t &g) {

  ssize_t i;
  a = c = t = g = 0;

  if (seq._ascii) {
    for (i = 0; i < seq.length; i++ ) {
      if (seq.seq[i] == 'a' || seq.seq[i] == 'A')
				a++;
      else if (seq.seq[i] == 'c' || seq.seq[i] == 'C')
				c++;
      else if (seq.seq[i] == 't' || seq.seq[i] == 'T')
				t++;
      else if (seq.seq[i] == 'g' || seq.seq[i] == 'G')
				g++;
    }
  }
  else {
    ssize_t nuc;
    for (i = 0; i < seq.length; i++ ) {
      nuc = DNASequence::UnmaskedValue(seq.seq[i]);
      if (nuc == 0)
				g++;
      else if (nuc == 1)
				a++;
      else if (nuc == 2)
				c++;
      else if (nuc == 3)
				t++;
    }
  }

}

ssize_t CountRepeatMasked(DNASequence &seq) {
  ssize_t i;
  ssize_t masked = 0;
  for (i = 0; i < seq.length; i++) 
		if (DNASequence::IsRepeatMaskedNuc(numeric_nuc_index[seq.seq[i]]))
			masked++;

  return masked;
}

void AppendReverseComplements(std::vector<SimpleSequence> &sequences) {
  ssize_t seq;
  ssize_t numSeq = sequences.size();
  SimpleSequence simpleSeq;
  simpleSeq.seq = NULL;
  sequences.resize(numSeq*2);
	// First stratify the sequences
	for (seq = numSeq-1; seq > 0; seq--) {
		sequences[seq*2].seq = sequences[seq].seq;
		sequences[seq*2].length = sequences[seq].length;
	}
	// fill in the reverse compliments between
  for (seq = 0; seq < numSeq; seq++ ) {
    MakeRC((char*) sequences[seq*2].seq, sequences[seq*2].length, sequences[seq*2+1].seq);
    sequences[seq*2+1].length = sequences[seq*2].length;
  }
}

void PrintSimpleSequences(std::vector<SimpleSequence> &sequences, std::string &outFileName) {
  std::ofstream out;
  openck(outFileName, out, std::ios::out);
  ssize_t i;
  DNASequence tmpSeq;
  tmpSeq._ascii = 1;
  for (i = 0; i < sequences.size(); i++) {
    tmpSeq.seq = sequences[i].seq;
    tmpSeq.length = sequences[i].length;
    tmpSeq.PrintSeq(out);
    *tmpSeq.titlestream << i;
    out << std::endl;
  }
  out.close();
}

ssize_t IsUnmasked(unsigned char ch) {
	if (std::isalpha((char)ch)) {
		return (ch == 'A' || ch == 'C' || ch == 'T' || ch == 'G');
	}
	else 
		return (ch < 4);
}
