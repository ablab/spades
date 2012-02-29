///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#ifndef FASTANAMEPARSER_H
#define FASTANAMEPARSER_H

#include "String.h"

class FastaNameParser {
 public:
  virtual ~FastaNameParser() {}
  virtual void extractNameFromBuffer( char* buffer, String& name ) const = 0;
};

class FullNameParser : public FastaNameParser {
 public:
  virtual ~FullNameParser() {}
  virtual void extractNameFromBuffer( char* buffer, String& name ) const;
};

class FirstWordParser : public FastaNameParser {
 public:
  virtual ~FirstWordParser() {}
  virtual void extractNameFromBuffer( char* buffer, String& name ) const;
};

class LastWordParser : public FastaNameParser {
 public:
  virtual ~LastWordParser() {}
  virtual void extractNameFromBuffer( char* buffer, String& name ) const;
};

class Riken_cDNA_Parser : public FastaNameParser {
 public:
  virtual ~Riken_cDNA_Parser() {}
  void extractNameFromBuffer( char* buffer, String& name ) const;
};

class TruncatedLastWordParser : public LastWordParser {
 public:
  virtual ~TruncatedLastWordParser() {}
  void extractNameFromBuffer( char* buffer, String& name ) const;
};


#endif

