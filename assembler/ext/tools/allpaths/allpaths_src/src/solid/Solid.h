/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef SOLID_H
#define SOLID_H

#include "Basevector.h"
#include "CoreTools.h"
#include "BasevectorTools.h"
#include "lookup/LookAlign.h"
#include "FastIfstream.h"
#include "TokenizeString.h"
#include <ext/hash_map>

using __gnu_cxx::hash_map;
using __gnu_cxx::hash;
struct eqstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return strcmp(s1, s2) == 0;
  }
};

typedef hash_map<const char*, unsigned int, hash<const char *>, eqstr> name2id_map_t;



template <typename T1, typename T2, typename T3>
struct trio {
  T1 first;
  T2 second;
  T3 third;

  trio(T1 v1, T2 v2, T3 v3) : first(v1), second(v2), third(v3) {}
  trio(const trio<T1,T2,T3> &t) : first(t.first), second(t.second), third(t.third) {}
  const trio<T1,T2,T3> & operator=(const trio<T1,T2,T3> &t) {
    first=t.first;
    second=t.second;
    third=t.third;
    return *this;
  }
};


template <typename T1, typename T2, typename T3>
inline trio<T1,T2,T3> make_trio(T1 v1, T2 v2, T3 v3) {
  return trio<T1,T2,T3>(v1,v2,v3);
}

typedef unsigned char SOLiD_color_t;
enum SOLiD_colors {
  SOLiD_BLUE = 0,
  SOLiD_GREEN = 1,
  SOLiD_ORANGE = 2,
  SOLiD_RED = 3
};

/// this table should be used to convert a pair
/// of adjacent bases to a SOLiD color (e.g.
/// AC converts to SOLiD_Bases2Color[BASE_A][BASE_C] = 
/// SOLiD_Bases2Color[0][1] = SOLiD_GREEN)
const unsigned char SOLiD_Bases2Color[4][4] = {
  { SOLiD_BLUE, SOLiD_GREEN, SOLiD_ORANGE, SOLiD_RED },
  { SOLiD_GREEN, SOLiD_BLUE, SOLiD_RED, SOLiD_ORANGE },
  { SOLiD_ORANGE, SOLiD_RED, SOLiD_BLUE, SOLiD_GREEN },
  { SOLiD_RED, SOLiD_ORANGE, SOLiD_GREEN, SOLiD_BLUE }
};


/// this table should be used to convert a base and
/// a color into the second base in the pair (conversion
/// into a color space is performed as <base1, base2> -> color,
/// so that to define the inverse conversion one has to know
/// the color and the first base in the pair: <base1,color> -> base2).
/// Usage: SOLiD_Color2Base[first_base][color] gives the base
/// that follows first_base in the base space sequence.
const unsigned char SOLiD_Color2Base[4][4] = {
  { BASE_A, BASE_C, BASE_G, BASE_T },
  { BASE_C, BASE_A, BASE_T, BASE_G },
  { BASE_G, BASE_T, BASE_A, BASE_C },
  { BASE_T, BASE_G, BASE_C, BASE_A }
};

/// Converts two bases into color representation
inline SOLiD_color_t SOLiD_toColor(base_t base1, base_t base2) {
  return SOLiD_Bases2Color[base1][base2];
}


/// Converts base and color into next base
inline base_t SOLiD_toBase(base_t base, SOLiD_color_t color) {
  return SOLiD_Color2Base[base][color];
}


/// \brief Converts vector of bases (sequence in base space)
/// into colorspace representation: first base, copied from
/// the original sequence, followed by sequence of SOLiD colors
/// representing pairs of adjacent bases.
///
/// @param bases [in] input sequence in base space
/// @param colors [out] colorspace representation of the input sequence

void ColorEncode( const basevector& bases, basevector& colors );


/// \brief Converts a sequence formed from <primer_base> followed
/// by <bases> 
/// into full colorspace representation: first base (i.e. <primer_base>) 
/// followed by sequence of SOLiD colors
/// representing pairs of adjacent bases.
///
/// @param primer_base full sequence to be converted is taken as <primer_base> followed by <bases>
/// @param bases [in] input sequence in base space
/// @param colors [out] full colorspace representation of the full input sequence (<primer_base>+<bases>)

void ColorEncode( base_t primer_base, const basevector& bases, basevector& colors );


/// \brief Converts a collection of base space sequences into
/// a collection of colorspace representations of these sequences.
///
/// Colorspace representation of the sequence: first base of the 
/// sequence followed by colors, each representing a pair of adjacent
/// bases in the original sequence.
void ColorEncode( const vecbasevector& bases, vecbasevector& colors );


/// \brief Converts a sequence from color space representation into
/// a base space representation. Colorspace representation is assumed
/// to have a base as its first element, followed by SOLiD colors,
/// representing adjacent base pairs. The base space representation will
/// start, accordingly, from the same base, followed by bases computed
/// by "unwinding" each preceding base and color back into a next base.
///
/// @param colors [in] input sequence in color space
/// @param bases [out] base space representation of the input sequence

void ColorDecode( const basevector & colors, basevector & bases );


/// \brief Converts a collection of color space sequences into
/// a collection of base space representations of these sequences.
///
/// Colorspace representation of the sequence must contain first base 
/// of the sequence followed by colors, each representing a pair of 
/// adjacent bases.

void ColorDecode( const vecbasevector& colors, vecbasevector& bases );

// wrapper class; currently used mainly to alleviate dispatching
// the methods capable of decorating their output with html bells and whistles
struct html_ostream {
  html_ostream(ostream & out) : _of(0), _o(out) {};
  html_ostream(const char * fname) :
    _of( new ofstream(fname) ),
    _o(*_of)
  {}

  ~html_ostream() { 
    if ( _of ) _of -> close();
    delete _of;
  }

  template <class T>
  friend html_ostream & operator << ( html_ostream & h, const T value ) {
    h._o << value;
    return h;
  }
private:
  ofstream *_of;
  ostream & _o;
};

struct tty_ostream  : public ostream {
  
};


/// Prints colorvector \c b into the output stream \c out adding VT220 
/// escape sequences to turn on colors. NOTE: as there is currently no difference
/// between internal representations of "true" basevectors and colorvectors,
/// this method \em assumes that the passed vector contains colors and prints
/// the data accordingly (as 0,1,2,3 rather than A,C,G,T). This method \em does
/// \em not expect the first element of the colorvector to be a base (followed by 
/// colors). If the colorvector used does contain a "full" colorspace sequence
/// (starting with a base), the first base should be removed manually (or 
/// argument \c start=1 should be used). This method does \em not print
/// newline character after the vector (think 'out << b;' - no newline).

void PrintAsColorsWithEsc(ostream & out, 
			  const basevector & b, 
			  unsigned int start = 0 , 
		       int length = -1 );

/// Prints colorvector \c b into the output stream \c out adding HTML 
/// tags to turn on colors. NOTE: as there is currently no difference
/// between internal representations of "true" basevectors and colorvectors,
/// this method \em assumes that the passed vector contains colors and prints
/// the data accordingly (as 0,1,2,3 rather than A,C,G,T). This method \em does
/// \em not expect the first element of the colorvector to be a base (followed by 
/// colors). If the colorvector used does contain a "full" colorspace sequence
/// (starting with a base), the first base should be removed manually (or 
/// argument \c start=1 should be used). This method does \em not print
/// newline/<br> character after the vector (think 'out << b;' - no newline).

void PrintAsColorsHTML(ostream & out, 
		       const basevector & b, 
			  unsigned int start = 0 , 
		       int length = -1 );


/// Prints colorvector \c b into the output stream. 
/// NOTE: as there is currently no difference
/// between internal representations of "true" basevectors and colorvectors,
/// this method \em assumes that the passed vector contains colors and prints
/// the data accordingly (as 0,1,2,3 rather than A,C,G,T). This method \em does
/// \em not expect the first element of the colorvector to be a base (followed by 
/// colors). If the colorvector used does contain a "full" colorspace sequence
/// (starting with a base), the first base should be removed manually (or 
/// argument \c start=1 should be used). This method does \em not print
/// newline character after the vector (think 'out << b;' - no newline).

void PrintAsColors(ostream & out, const basevector & b, unsigned int start = 0 , int length = -1 );


/// Prints reverse of the colorvector \c b (i.e. colorvector that would
/// result from the reverse complement of the original basespace sequence)
/// into the output stream. 
/// NOTE: as there is currently no difference
/// between internal representations of "true" basevectors and colorvectors,
/// this method \em assumes that the passed vector contains colors and prints
/// the data accordingly (as 0,1,2,3 rather than A,C,G,T). This method \em does
/// \em not expect the first element of the colorvector to be a base (followed by 
/// colors). If the colorvector used does contain a "full" colorspace sequence
/// (starting with a base), the first base should be removed manually (or \c length
/// argument can be set appropriately). \c start position is counted from the 
/// \em end of the vector \b (\em i.e. from the \em beginning of the reverse
/// colorvector to be printed). This method does \em not print
/// newline character after the vector (think: out << b; - no newline).
void PrintReverseAsColors(ostream & out, const basevector & b, unsigned int start = 0 , int length = -1 );



/// Prints visual alignment of colorvector b. Alignment is specified by the
/// look_align argument and the basespace reference is passed in ref argument.
/// b is assumed to be in pure colorspace representation; if the first base is known
/// it should be passed separately in <first_base> argument  
void PrintVisualColorAlign(ostream & out, const basevector & b, 
			    const vecbasevector & ref,
			    const look_align & la, 
			   base_t first_base = 255);

/// Builds full colorspace representation (with the leading primer base \c first_base) of the target sequence
/// pointed to by the alignment \c la on the basespace reference \c ref (the sequence is also properly reverse 
/// complemented is the alignment is reverse). In other words, if you have a (full) colorspace read and its alignment,
/// this method will extract the reference sequence (in proper orientation) the read's DNA molecule purportedly
/// represents. The resulting colorspace sequence is
/// written into \c result (previous content, if any, is overwritten and the vector is resized)
void FullColorspaceTargetSeq(const look_align &la, const vecbasevector &ref, base_t first_base, basevector &result);


/// Applies simple one pass correction to the full colorspace reference \c read aligned as specified by \c la
/// to the reference \c ref [note: the reference must be provided in the original basespace representation]. The
/// error-corrected full colorspace sequence of the read will be written into \c result (previous content, if any,
/// is destroyed, and the vector is resized).
///
/// One-pass correction strategy is as follows: 1) every single mismatch flanked by matches is corrected;
/// 2) mismatch followed by a discordant mismatch is immediately corrected; 3) mismatch (1) followed by a concordant
/// mismatch (2) is assumed to be a SNP and never corrected, the position (2) is treated as match for the purposes
/// of correcting subsequent colors. Note that when this strategy is used, then, for instance, a sequence of
/// 3 mismatches (1)(2)(3) such that both (1),(2) and (2),(3) are concordant will be always called as a SNP at (1)(2)
/// (read colors (1)(2) will not be corrected) and a mismatch at (3) (i.e. color at (3) will get corrected).
///
/// This method returns a trio of counts: \c first is the number of (corrected) mismatches (those that did not end up in SNPs), 
/// \c second is the number of discordant pairs (\em adjacent pairs are counted, i.e. *** is 1 pair and **** is 2 pairs),
/// and \c third is the number of SNPs (concordant mismatches in two adjacent positions (pairs)) chosen in one pass as 
/// described above, so that total number of mismatching colors is \c first + 2*\c third.
trio<unsigned int, unsigned int,unsigned int> CorrectColorRead( const basevector & read, const vecbasevector & ref, 
						   const look_align & la, basevector & result );


/// Load generic SOLiD (abi) aligns into the collector; <nmap> provides the mapping
/// from read names to read ids. Only aligns for the reads that have their names
/// stored in the map will be extracted; query_id for the aligns of there reads will be set to the mapped values.
template <class Collector>
void LoadSOLiDAligns(const String & fname, 
		     const name2id_map_t & nmap,
		     Collector & aligns) {

  fast_ifstream in(fname);

  vec<char> align_separator;
  align_separator.push_back(',');
  const String tig_pos_sep('_');
  const String pos_err_sep('.');

  look_align la;
  la.nhits=la.indels=0;
  la.a.SetNblocks(1);
  la.a.SetGap(0,0);
  la.SetStartOnQuery(0);

  String header_line;
  String seq_line;
  //  int read_id = -1;
  int cnt = 0;

  name2id_map_t::const_iterator map_end = nmap.end();

  while  (1) {
    getline(in, header_line);
    if ( in.fail() ) break;
    getline(in, seq_line);
    if ( in.fail() ) break;

    if ( header_line[0] != '>' ) {
      cout << "Could not find alignment header where expected"<< endl;
      exit(1);
    }
    
    //    read_id++;
    //    if ( low_id > 0 && read_id < low_id ) continue;
    //    if ( high_id > 0 && read_id >= high_id ) break;    
    
    vec<String> tokens;
    Tokenize(header_line,align_separator,tokens);
    //    cout << tokens[0] << ':' << tokens.size() << endl;
    if ( tokens.size() == 1 ) continue; //no aligns
    //    if ( tokens[0] == ">198_449_1400_R3") {
    //      cout << endl << "Found what we know as 12970980" << endl;
    //          PRINT( read_id );
    //          cout << header_line << endl;
    //     }
    
    name2id_map_t::const_iterator it = nmap.find( tokens[0].c_str()+1 );
    if ( it == map_end ) continue; // the read name is not among the set of acceptable names

    la.query_id = it->second;
    seq_line.Trim(" \t");
    la.query_length = seq_line.size() - 2;
    la.a.SetLength(0,la.query_length);

    for ( unsigned int i = 1 ; i < tokens.size() ; i++ ) {
      la.target_id = tokens[i].Before(tig_pos_sep).Int();
      int pos = tokens[i].Between(tig_pos_sep,pos_err_sep).Int();
      if ( pos < 0 ) {
	  la.rc1 = true;
	  pos = (-pos);
	  //	  la.a.Setpos2(pos-la.query_length);
	  la.SetStartOnTarget(pos-la.query_length);
      } else {
	  la.rc1 = false;
	  //	  la.a.Setpos2(pos);
	  la.SetStartOnTarget(pos);
      }
      la.mutations = tokens[i].After(pos_err_sep).Int();
      //      if ( read_id == 12970980 ) {
      //	cout << endl << "found putative read 12970980 (according to solid):" << endl;
      //      	la.PrintReadableBrief(cout);
      //      	cout << header_line << endl;
      //      }
      aligns.Insert(la);
    }

  }

  aligns.Consolidate();
}



#endif
