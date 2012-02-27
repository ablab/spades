/***************************************************************************
 * Title:          DNASequence.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/26/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "DNASequence.h"
#include "SeqUtils.h"
#include "compatibility.h"

char DNASequence::_nucs[26] = { 1,12, 2,14,100,100, 0,11,100,100,10,100, 9, 4,100,100,100, 5, 8, 3,100,13, 7,15, 6,100};

/*
char DNASequence::_indexToNuc[16] = {'G', 'A', 'C', 'T', 
				     'N', 'R', 'Y', 'W', 
				     'S', 'M', 'K', 'H'}
*/				     


ssize_t DNASequence::FORWARD_STRAND = 0;
ssize_t DNASequence::REVERSE_STRAND = 1;

char DNASequence::_nucToIndex[256] = {16,16,16,16,16,16,16,16,16,16,  // 0
				      16,16,16,16,16,16,16,16,16,16,  // 10
				      16,16,16,16,16,16,16,16,16,16,  // 20
				      16,16,16,16,16,16,16,16,16,16,  // 30
				      16,16,16,16,16,16,16,16,16,16,  // 40
				      16,16,16,16,16,16,16,16,16,16,  // 50
				      16,16,16,16,16,0,12,1,14,16,  // 60
				      16,2,11,16,16,116,16,9,4,16,  // 70
				      16,16,5,8,3,16,13,7,15,6,  // 80
				      16,16,16,16,16,16,16,252,12,253,  // 90  251..255 are the masked characters.
				      14,16,16,254,11,8,16,116,16,9,  // 100
				      4,16,16,16,5,16,255,16,13,7,  // 110
				      15,6,16,16,16,16,16,16,16,16,  // 120
				      16,16,16,16,16,16,16,16,16,16,  // 130
				      16,16,16,16,16,16,16,16,16,16,  // 140
				      16,16,16,16,16,16,16,16,16,16,  // 150
				      16,16,16,16,16,16,16,16,16,16,  // 160
				      16,16,16,16,16,16,16,16,16,16,  // 170
				      16,16,16,16,16,16,16,16,16,16,  // 180
				      16,16,16,16,16,16,16,16,16,16,  // 190
				      16,16,16,16,16,16,16,16,16,16,  // 200
				      16,16,16,16,16,16,16,16,16,16,  // 210
				      16,16,16,16,16,16,16,16,16,16,  // 220
				      16,16,16,16,16,16,16,16,16,16,  // 230
				      16,16,16,16,16,16,16,16,16,16,  // 240
				      16,16,16,16,16,16};         // 250

char DNASequence::_indexToNuc[256] = {'A', 'C', 'G', 'T', 'N', 'R', 'Y', 'W', 'S', 'M',   // 9
				      'K', 'H', 'B', 'V', 'D', 'X', '\0','\0','\0','\0',  // 19
				      '\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',  // 29 
				      '\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',  // 39 
				      '\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',  // 49 
				      '\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',  // 59 
				      '\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',  // 69 
				      '\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',  // 79 
				      '\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',  // 89 
				      '\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',  // 99
				      '\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',  // 109
				      '\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',  // 119
				      '\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',  // 129
				      '\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',  // 139
				      '\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',  // 149
				      '\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',  // 159
				      '\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',  // 169
				      '\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',  // 179
				      '\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',  // 189
				      '\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',  // 199
				      '\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',  // 209
				      '\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',  // 219
				      '\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',  // 229
				      '\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',  // 239
				      '\0','\0','\0','\0', 'g', 'a', 'c', 't', 'g', 'a',  // 249
				      'c', 't','g','a','c','t'};                          // 255


char DNASequence::CharToNumeric(char nuc) {
  ssize_t index = (ssize_t) nuc;
  //UNUSED// ssize_t translation;
  // masked sequence
  if ((_masked & REPEAT_MASKED)){
    return _nucToIndex[index];
  }
  else {
    return UnmaskedValue(_nucToIndex[index]);
  }
}

void DNASequence::MaskPosition(ssize_t pos) {
  assert(seq);
  assert(pos < length);
  if (_ascii == 0) {
    seq[pos] = -abs(seq[pos]);
  }
  else {
    seq[pos] = tolower(seq[pos]);
  }
}
  
const unsigned char & DNASequence::operator[](const ssize_t &pos) {
  assert (pos >= 0 && pos < length);
  return seq[pos];
}

char DNASequence::Get2BitValue(char numericNuc) {
  assert(numericNuc >= -12 && numericNuc <= 3);
  return UnmaskedValue(numericNuc);
}

char DNASequence::GetNucChar(char numericNuc) {
  ssize_t index = (unsigned char) numericNuc;
  return _indexToNuc[index];
}
  
bool DNASequence::ValidNuc(char nuc, ssize_t maskType) {
  // Non-standard base here.  Anything that is invalid here is invalid for 
  // all types of masking.
  if (UnmaskedValue(nuc) >= 4) 
    return 0;

  if (maskType & REPEAT_MASKED) {
    if (IsRepeatMaskedNuc(nuc)) 
      return 0;
  }
  
  if (maskType & UNIQUE_MASKED) {
    if (IsUniqueMaskedNuc(nuc)) 
      return 0;
  }
  
  if (maskType & COMMON_MASKED) {
    if (IsCommonMaskedNuc(nuc)) 
      return 0;
  }
  // Ok to here
  return 1;
}

bool DNASequence::ValidSequence(ssize_t pos, ssize_t searchLen, ssize_t allowMask) {
  ssize_t i;
  for (i = pos; i < pos + searchLen; i++ ) {
    if (! ValidNuc(seq[i], allowMask))
      return 0;
  }
  return 1;
}

void DNASequence::StoreAsMasked(ssize_t masked) {
  _masked = masked;
}
  
std::ostream & DNASequence::PrintASCIISeq(std::ostream &out, ssize_t lineLength) {

	ssize_t curPos = 0;
	ssize_t curLineLength;
	char newline = '\n';
	while (curPos < length) {
		if (length - curPos > lineLength)
			curLineLength = lineLength;
		else 
			curLineLength = length - curPos;
		out.write((char*) &seq[curPos], curLineLength);
		curPos += curLineLength;
		if (curPos < length)
			out.write(&newline, 1);

	}
	return out;
}

std::ostream & DNASequence::PrintSeq(std::ostream &out, ssize_t lineLength) {
	//  long i, j;
	//  i = 0;
	_SZT_ j;
  if (namestr != "" && lineLength > 0)
    out << ">" << namestr << std::endl;
  else if (titlestream->str() != "") {
    out << ">" << titlestream->str() << std::endl;
    namestr = "";
    titlestream->str(namestr);
  }
	if (_ascii) {
		// Fast sequence writer.
		PrintASCIISeq(out, lineLength);
		return out;
	}
  for (j = 0; j < length ; j++) {
    if (!_ascii && (seq[j] < 243 and seq[j] > 16)) {
      std::cout << "error in seqn  " << seq[j] << " at position: " << j << std::endl;
      exit(0);
    }
    if (_ascii)
      out << seq[j];
    else {
      out << GetNucChar(seq[j]);
    }
    if (lineLength > 0 && j % lineLength == (lineLength-1) && j != length) 
      out << std::endl;
  }
  return out;
}

std::ostream &DNASequence::PrintlnSeq(std::ostream &out, ssize_t lineLength) {
	PrintSeq(out, lineLength);
	out << std::endl;
	return out;
}

void DNASequence::MergeSequence(DNASequence &other, ssize_t mergeMasked) {

  assert("must merge same lengths" && length == other.length);

  ssize_t pos;
  for (pos = 0; pos < length; pos++) {
    if (mergeMasked) {
      // masked nucleotides in other take precedence
      // over nucleotides in this sequence
      if (other.IsPosRepeatMasked(pos)) {
	MarkPosRepeatMasked(pos);
      }
      if (other.IsPosUniqueMasked(pos)) {
	MarkPosNotUnique(pos);
      }
      if (other.IsPosCommonMasked(pos)) {
	MarkPosNotCommon(pos);
      }
    }
    else {
      // unmasked nucleotides in 'other' take precedence 
      // over nucleotides in this sequence
      if (! other.IsPosMasked(pos) ) {
	UnmaskPos(pos);
      }
    }
  }
}


void DNASequence::MergeSequence(char* src, ssize_t srcLength) {

  assert("must merge same lengths" && length == srcLength);
  DNASequence other;
  other.InitializeNoCopy(src, srcLength);
  ssize_t pos;
  for (pos = 0; pos < length; pos++) {
    if (other.IsPosRepeatMasked(pos)) {
      MarkPosRepeatMasked(pos);
    }
    if (other.IsPosUniqueMasked(pos)) {
      MarkPosNotUnique(pos);
    }
    if (other.IsPosCommonMasked(pos)) {
      MarkPosNotCommon(pos);
    }
  }
}

void DNASequence::HardMask() {
  // This only works on unconverted sequences
  ssize_t i;
  for (i = 0; i < length; i++) {
    char nuc;
    nuc = seq[i];
    if (nuc >= 'a' and nuc <= 'z')
      seq[i]= 'N';
  }
}

void DNASequence::RemoveRepeats() {
	ssize_t i;
	for (i = 0; i < length; i++) {
		seq[i] = unmasked_nuc[seq[i]];
	}
}
