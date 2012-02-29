///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef EFASTA_TOOLS_H
#define EFASTA_TOOLS_H

#include "Basevector.h"
#include "Bitvector.h"
#include "CoreTools.h"
#include "Fastavector.h"
#include "Superb.h"

struct Ambiguity
{
  size_t start;
  int size;
  String replace;

  String to_annotation() const 
  {
    String s = "";
    if (size == 0) {
      s += (ToString(start) + "^\t" +
            ToString(start + 1) + "\tvariation\t\t\n");
    }
    else {
      s += (ToString(start + 1) + "\t" +
            ToString(start + size) + "\tvariation\t\t\n");
    }
       
    s += "\t\t\treplace\t" + ((replace == "") ? "\"\"" : replace) + "\n";
    
    s += "\t\t\tnote\tassembly ambiguity\n"; 
    return s;
  }
};

// Class efasta represents a single efasta record, exclusive of the header line.

class efasta : public String {

     public:

     efasta( ) { }
     efasta( const String& s ) : String(s) { }
     efasta( const fastavector& v );  
     efasta( const vec<basevector>& x );

     // Length1 returns the length in bases of the record, assuming that the
     // first choice is made for all brackets.

     int Length1( Bool count_Ns = false ) const;

     // At what index in the efasta does the nth base occur, assuming first choices?
     // If it's in an alternative, use its beginning, because we may be doing this in
     // order to chop the sequence here (could be handled better). --bruce

     int Index1(const int n) const;

     // Index1Alt: similar to Index1, but does a direct translation.  For example,
     // 0123456
     // C{A,T}G
     // C A   G
     // 0 1   2
     // 0 maps to 0, 1 maps to 2, 2 maps to 6.

     int Index1Alt(const int n) const;

     // MinLength returns the shortes length in bases of the record
     
     int MinLength( ) const;

     
     // MaxLength returns the largest length in bases of the record
     
     int MaxLength( ) const;

     // Ambiguities returns the number of "ambiguities" in the record, which we
     // define to be the sum of all expressions {x1,...,xn}, of n-1.

     int Ambiguities( ) const;

     // AmbCount returns the number of "ambiguous bases" in the record, which we
     // define to be the sum over all expressions {x1,...,xn}, of the maximum
     // of the lengths of the xi.

     int AmbCount( ) const;
     int AmbCount(int& snp_count, int&indel_count ) const;

     // Indicate if the index position n is within brackets (in the ambiguous region) 

     Bool IndexInAmbiguity( const longlong pos ) const;

     // FlattenTo.  Make the first choice when a bracket expression is encountered.
     // When converting to a basevector, converts Ns to As.  The form
     // FlattenTo(fastavector) converts expressions like {A,C} into ambiguous bases.

     void FlattenTo( basevector& b ) const;
     void FlattenTo( basevector& b, bitvector& gaps ) const;
     void FlattenTo( fastavector& v ) const;
     void FlattenTo( fastavector& v, vec<Ambiguity> & va, 
          const Bool ambiguous_base_codes = True) const;

     // FlattenMaxTo. When a bracket expression is encountered, make the longest 
     // choice but convert each character to 'N'.

     void FlattenNMaxTo( fastavector& v ) const;

     // ExpandTo.  Convert to a list of fastavectors.  This will find ambiguous
     // base codes.   If max_count is specified and the number of fastavectors would 
     // exceed it, leave v empty and return False.

     Bool ExpandTo( vec<fastavector>& v, const int max_count = -1 ) const;
     Bool ExpandTo( vec<basevector>& v, const int max_count = -1 ) const;

     // MakeGraph.  Convert to a graph.  Note that in the basevector case, this 
     // will convert Ns to As.

     void MakeGraph( vec< vec<basevector> >& G ) const;
     void MakeGraph( vec< vec<fastavector> >& G ) const;

     void MakeFromGraph( const vec< vec<basevector> >& G );

     // GetKmers.  Find some (or all) of the kmers, appending.  If max_per is 
     // specified, the number of kmers starting at a given position is capped at 
     // that value.  The exact set of kmers that are returned in that case is 
     // determined by what the code does.  Note that not setting max_per is 
     // dangerous.  This code will return duplicates in certain instances.  For 
     // example, CT{T,}A will yield the 2-mer TA twice.  It is not exactly clear 
     // exactly how duplicates should be removed, as for example CTTTA legitimately 
     // has the 2-mer TT twice.

     void GetKmers( const int K, vecbasevector& kmers, 
          const int max_per = -1 ) const;

      // Prints in a fasta format: "><string_id>\n" followed by the full base
      // sequence stored in the basevector; breaks
      // long sequences nicely into 80-character lines

      void Print( ostream& out, const String& id = "" ) const;
 
      // Write a efasta-format file which contains the bases in the
      // input fasta scaffolded together (with gaps, etc.) as defined
      // by the scaffolds.  If supplied, rc defines which contigs are
      // reverse-complement.

      friend void WriteScaffoldedEFasta( const String &out_file,
					const vec<efasta> &fasta,
					const vec<superb> &scaffolds,
					const Bool ncbi_format = False);
    
};

// ExpandAmbCode.  Expand ambiguous base codes.  For example, M is converted to
// {A,C}, and N is converted to {A,C,G,T}.  This also converts n to N.

String ExpandAmbCode( char x );
String ExpandAmbCode( const String& x );

// ValidateEfastaRecord.  Test a group of lines to see if they represent a valid 
// efasta record, exclusive of the header line.  This doesn't check for duplicates 
// in choose expressions.

void ValidateEfastaRecord( const vec<String>& lines, const String& msg = "" );
void ValidateEfastaRecord( const String& line, const String& msg = "" );

// LoadEfastaFlat.  Load an efasta file, making the first choice when a
// bracket expression is encountered.  This converts Ns to As.

void LoadEfastaFlat( const String& fn, vecbasevector& bases );

// LoadEfastaFlatGaps.  Load an efasta file, taking the first choice
// when a bracket expression is encountered, and returning the location of gaps.
// This is parallel to LoadEfastaFlat.

void LoadEfastaFlatGaps( const String& fn, vecbitvector& gaps );

// LoadEfastaIntoStrings.  Load an efasta file, yielding a vector of strings, one 
// for each record.  Newlines are removed.  Header lines are discarded.

void LoadEfastaIntoStrings( const String& fn, vec<efasta>& x );

// SplitEfastaIntoContigs.  Given an efasta file that has been loaded into a
// vector of efastas, split it into contigs, yielding as output a new vector of
// efastas (one for each contig) and a vec<superb> that describes the scaffold 
// structure.  Contig lengths are based on the first choice in each bracket.
// Gap deviations are set to zero.

void SplitEfastaIntoContigs( const vec<efasta>& scaffolds,
     vec<efasta>& contigs, vec<superb>& scaffold_structure );

// AllPlus.  Return the concatenation of some lines, including newlines.

String AllPlus( const vec<String>& lines );

// ConvertToLengths.  Replace each instance of {x1,...,xn} by
// {s1,...,sn} where the si are the lengths of the xi, sorted.

String ReplaceByLengths( const efasta& e );

typedef efasta EFastaVec;

// Given a collection of sequences to be turned into efasta, find the number of 
// bases shared by the patches have on the left.  Then, excluding the shared bases 
// that have already been declared on the left, determine the number of shared bases 
// that the patches have on the right.

void GetShares( const vec<basevector>& patches, int& left_share, int& right_share );


#endif
