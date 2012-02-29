///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#ifndef FASTAVERIFIER_H
#define FASTAVERIFIER_H

// FastaVerifier understands the fasta format and can confirm whether
// a given string is a valid line in that format.

class FastaVerifier {

 public:
  virtual ~FastaVerifier() {}

  bool verifyLine( const char* line );
  bool verifyRestOfLine( const char* line );
  
 protected:
  virtual bool verifyChar( const char c ) = 0;
  
  bool ignore_rest_of_line_;
};

class FastaNullVerifier: public FastaVerifier {
 
  private:
    bool verifyChar( const char c ) { return true; }
};

class FastaSequenceVerifier: public FastaVerifier {
 
 private:
  bool verifyChar( const char c );
};

class FastaQualityVerifier: public FastaVerifier {

 private:
  bool verifyChar( const char c );
};

#endif
