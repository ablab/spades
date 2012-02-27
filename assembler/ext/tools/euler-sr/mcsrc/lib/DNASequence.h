/***************************************************************************
 * Title:          DNASequence.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/28/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef _DNA_SEQUENCE 
#define _DNA_SEQUENCE

#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <ostream>
#include <sstream>

#include "SimpleSequence.h"
#include "utils.h"

#define IGNORE_MASKED 0
#define REPEAT_MASKED 1
#define UNIQUE_MASKED 2
#define COMMON_MASKED 4

#ifndef getmin
#define getmin(a,b) (((a) < (b)) ? (a) : (b))
#endif

#ifndef getmax
#define getmax(a,b) (((a) > (b)) ? (a) : (b))
#endif
/*
class CopyableStream {
public:
  std::stringstream strm;
  CopyableStream & operator=(const CopyableStream &copy) {
    if (this != &copy) {
      std::string strcopy;
      strcopy = copy.str();
      this->strm.str(strcopy);
    }
    return *this;
  }

};
*/

class NamedSequence : public SimpleSequence {
 public:
	std::string namestr;
};

class DNASequence : public NamedSequence {
public:
  ssize_t  strand;
  ssize_t  _masked;
  ssize_t _ascii;
  static char nucChar[256];
  static char _nucs[26];
  static char _indexToNuc[256];
  static char _nucToIndex[256];
  ssize_t startEnumeration;
  ssize_t startPosition;
  std::stringstream *titlestream;
  std::string namestr; // TODO: 'namestr' is also ~15 lines above in base class NamedSequence; do we need it twice or should we get rid of this one? 
  static ssize_t FORWARD_STRAND;
  static ssize_t REVERSE_STRAND;

  DNASequence() {
    seq = NULL;
    length = 0;
    _ascii = 1;
    _masked = 1;
    namestr = "";
    titlestream = new std::stringstream;
    titlestream->str(namestr);
    strand = FORWARD_STRAND;
  }
  
	// TODO: Fix sparc compiler warning
	// "warning: base class `class NamedSequence' should be explicitly initialized in the copy constructor"
	// First deal with double definition of 'namestr', then init it correctly
  DNASequence(const DNASequence &copy) {
    titlestream = new std::stringstream;
    *this = copy;
  }

  DNASequence(ssize_t newLength) {
    namestr = "";
    titlestream->str(namestr);
    seq  = new unsigned char[newLength];
    length = newLength;
    memset(seq, 0, newLength);
    _ascii = 1;
    startPosition = 0;
    startEnumeration = 0;
    strand = FORWARD_STRAND;
  }
  
  DNASequence& Reference(std::string &str) {
    seq = (unsigned char*) str.c_str();
    length = str.size();
    return *this;
  }

  DNASequence& CopyString(std::string &str) {
    Reset(str.size());
    memcpy(seq, (unsigned char*) str.c_str(), str.size());
    length = str.size();
		return *this;
  }

  DNASequence &operator=(const SimpleSequence &src) {
		if (this != &src) {
			seq = src.seq;
			length = src.length;
		}
    return *this;
  }

  DNASequence &operator=(const DNASequence &src) {
		if (this == &src)
			return *this;

    if (seq != NULL)
      delete[] seq;
    seq = new unsigned char[src.length];
    memcpy(seq, src.seq, src.length);
    length = src.length;
    strand = src.strand;

    namestr          = src.namestr;
    titlestream->str(namestr);
    _ascii           = src._ascii;
    startEnumeration = src.startEnumeration;
    startPosition    = src.startPosition;
    return *this;
  }

  ~DNASequence() {
    delete titlestream;
  }

  void SetAscii() { _ascii = 1; }

  void InitializeCopy(char *src, ssize_t srcLength) {
    if (seq != NULL) 
      delete seq;
    
    seq = new unsigned char[srcLength];
    memcpy((void*)seq, (const void*)src, srcLength);
  }

  void ClearName() {
    namestr = "";
    titlestream->str(namestr);
  }

  void StoreName(char *srcName) {
    namestr = srcName;
  }
  void StoreName(std::string newname) {
    namestr = newname;
  }

  void InitializeNoCopy(char *src, ssize_t srcLength) {
    if (seq != NULL) 
      delete seq;

    seq = (unsigned char*) src;
    length = srcLength;
  }
  
  void Free() {
    if (seq != NULL) {
      delete [] seq;
      seq = NULL;
    }
  }
  void Reset(ssize_t newLength=0) {
    if (seq != NULL)
      delete[]  seq;
    if (newLength != 0) {
      seq = new unsigned char[newLength];
      memset(seq, 0, newLength);
    }
    else 
      seq = NULL;
    length = newLength;
  }
  void CopyDetails(DNASequence &seq) {
    _ascii = seq._ascii;
    _masked = seq._masked;
    namestr  = seq.namestr;
  }
  void Copy(DNASequence &src, ssize_t start = -1, ssize_t end = -1) {
    CopyDetails(src);
    if (start == -1) {
      Reset(src.length);
    }
    else {
      assert(end != -1);
      Reset(end - start);
      memcpy(this->seq, &src.seq[start], end - start);
    }
  }

  void MaskPosition(ssize_t pos);

  const unsigned char & operator[](const ssize_t &pos);

  // Functions for dealing with representations of nucleotides.
  static char Get2BitValue(char numericNuc);
  static char GetNucChar(char numericNuc);

  ssize_t GetPos(ssize_t pos, ssize_t lookupStrand) {
    if (lookupStrand == FORWARD_STRAND) {
      if (strand == FORWARD_STRAND) {
	return pos;
      }
      else {
	return (length - pos - 1);
      }
    }
    else {
      if (strand == FORWARD_STRAND) {
	return (length - pos - 1);
      }
      else {
	return pos;
      }
    }
  }

  char CharToNumeric(char nuc);
  // This will probably be deprecated
  bool ValidNuc(char nuc, ssize_t allowMask);

  // This will probably be deprecated as well.
  // Check to see if none of a sequence is masked.
  bool ValidSequence(ssize_t pos, ssize_t searchLen, ssize_t allowMask);

  // Pritnt a formated version of this sequence to stdout.
  // It includes masking information.
	std::ostream& PrintASCIISeq(std::ostream &out, ssize_t lineLength);
  std::ostream& PrintSeq(std::ostream &out, ssize_t length = 60);
	std::ostream& PrintlnSeq(std::ostream &out, ssize_t length = 60);
  // masked = 1 means that repmask masked nucleotides are stored as masked.
  void StoreAsMasked(ssize_t masked);

  // Convention for storing nucleotides:
  // The states of nucleotides are masked by different values.  They may 
  // be:
  //       Category                                       value (g,a,c,t)
  //  ------------------------------------------------------------
  //   1.  Normal                                         0,1,2,3 (g,a,c,t)
  //   2.  Repeat masked                                  -4,-3,-2,-1 (t,c,a,g)
  //   3.  Not unique in any number of sequences.         -8,-7,-6,-5 (t,c,a,g)
  //   4.  Not common between any number of sequences.    -12,-11,-10,-9 (t,c,a,g)
  //

	void RemoveRepeats();

  void MergeSequence(DNASequence &seq, ssize_t maskDown = 1);
  void MergeSequence(char* seq, ssize_t length);

  static char UnmaskedValue(char numericNuc) { 
    if (numericNuc < 0) 
      return (-numericNuc -1) % 4;
    else
      return numericNuc;
  }

  // Utility funcitons for computing states.
  static char ToRepeatMasked(char numericNuc) {return -UnmaskedValue(numericNuc) - 1;}
  static char ToNotUnique(char numericNuc)    {return -UnmaskedValue(numericNuc) - 5;}
  static char ToNotCommon(char numericNuc)    {return -UnmaskedValue(numericNuc) - 9;}

  // Check for general nucleotide mask
  // Return true if the nucleotide is masked at all
  static bool IsMasked(char numericNuc) { return (!(numericNuc >= 0 && numericNuc < 4)); }
  
  // Check for specific mask set on a nucleotide.
  static bool IsRepeatMaskedNuc(char numericNuc) {return (numericNuc <= -1 && numericNuc >= -4); }
  static bool IsUniqueMaskedNuc(char numericNuc) {return (numericNuc <= -5 && numericNuc >= -8); }
  static bool IsCommonMaskedNuc(char numericNuc) {return (numericNuc <= -9 && numericNuc >= -12); }

  ssize_t  IsACTG(ssize_t pos) {
    return UnmaskedValue(seq[pos]) < 4;
  }

  char MarkPosRepeatMasked(ssize_t pos) {
    assert(pos < length); 
    //    assert(UnmaskedValue(seq[pos]) < 4); 
    if (UnmaskedValue(seq[pos]) < 4)
      return seq[pos] = ToRepeatMasked(seq[pos]);
    else 
      return seq[pos];
  }
  char MarkPosNotUnique(ssize_t pos) {
    assert(pos < length); 
    //    assert(UnmaskedValue(seq[pos]) < 4);  
    if (UnmaskedValue(seq[pos]) < 4)
      return seq[pos] = ToNotUnique(seq[pos]);
    else 
      return seq[pos];
  }
  char MarkPosNotCommon(ssize_t pos){
    assert(pos < length); 
    if(UnmaskedValue(seq[pos]) < 4) 
      return seq[pos] = ToNotCommon(seq[pos]);
    else
      return seq[pos];
  }

  char UnmaskNoCk(ssize_t maskStart = 0, _SSZT_ maskLength=-1) {
    ssize_t i;
    if (maskLength == -1)
      maskLength = length;
    maskLength = std::min(length - maskStart, maskLength);
    for(i = maskStart; i < maskLength; i++) {
      UnmaskPosNoCk(i);
    }
    return 0;
  }
    
  char Unmask() {
    ssize_t i;
    for(i = 0; i < length; i++) {
      UnmaskPos(i);
    }
    return 0;
  }
  
  char UnmaskPos(ssize_t pos) {
    assert(pos < length);
    assert(UnmaskedValue(seq[pos]) < 4);
    return seq[pos] = UnmaskedValue(seq[pos]);
  }

  char UnmaskPosNoCk(ssize_t pos) {
    return seq[pos] = UnmaskedValue(seq[pos]);
  }
  // Check for specific mask set on a nucleotide.
  bool IsPosRepeatMasked(ssize_t pos) {assert(pos< length); return IsRepeatMaskedNuc(seq[pos]);}
  bool IsPosUniqueMasked(ssize_t pos){assert(pos< length); return IsUniqueMaskedNuc(seq[pos]);}
  bool IsPosCommonMasked(ssize_t pos){assert(pos< length); return IsCommonMaskedNuc(seq[pos]);}
  bool IsPosMasked(ssize_t pos) {assert(pos < length); return IsMasked(seq[pos]);}

  ssize_t CoordInRC(ssize_t pos) {
    return length - pos + 1;
  }
  void HardMask();
	void AddStrandKey(ssize_t strand) {
		namestr += " strand=" + NumToStr(strand);
	}
	void AddPosKey(ssize_t pos) {
		namestr += " pos=" + NumToStr(pos);
	}

};

typedef std::vector<DNASequence> DNASequenceList;
typedef std::vector<NamedSequence> NamedSequenceList;
#endif
