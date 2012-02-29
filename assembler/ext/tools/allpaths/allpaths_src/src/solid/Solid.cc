/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "CoreTools.h"
#include "solid/Solid.h"


/// \brief Converts vector of bases (sequence in base space)
/// into colorspace representation: first base, copied from
/// the original sequence, followed by sequence of SOLiD colors
/// representing pairs of adjacent bases.
///
/// @param bases [in] input sequence in base space
/// @param colors [out] colorspace representation of the input sequence

void ColorEncode( const basevector& bases, basevector& colors )
{   
  colors.resize(bases.size()); // reserve space
  if ( bases.size() == 0 ) return; // nothing to do
  basevector::const_iterator b_iter = bases.Begin();
  base_t base1 = *b_iter; ++b_iter; // get the first base
  // and store the first base into the colorspace representation:
  colors.Set(0, base1); 

  // now run through the remaining bases and convert each pair of
  // adjacent bases into a color (base1, base2) --> color
  for ( unsigned int j = 1; j < bases.size( ); j++ ) {
    base_t base2 = *b_iter; ++b_iter; // get next base
    // convert (base1, base2) into a color:
    colors.Set( j, SOLiD_toColor( base1, base2) );
    base1 = base2; // shift and prepare to retrieve the next base
  }
}

/// \brief Converts a sequence formed from <primer_base> followed
/// by <bases> 
/// into full colorspace representation: first base (i.e. <primer_base>) 
/// followed by sequence of SOLiD colors
/// representing pairs of adjacent bases.
///
/// @param primer_base full sequence to be converted is taken as <primer_base> followed by <bases>
/// @param bases [in] input sequence in base space
/// @param colors [out] full colorspace representation of the full input sequence (<primer_base>+<bases>)

void ColorEncode( base_t primer_base, const basevector& bases, basevector& colors )
{   
  colors.resize(bases.size()+1); // reserve space
  base_t base1 = primer_base;
  // store the first base into the colorspace representation:
  colors.Set(0, base1); 

  basevector::const_iterator b_iter = bases.Begin();
  
  // now run through the remaining bases and convert each pair of
  // adjacent bases into a color (base1, base2) --> color
  for ( unsigned int j = 1; j <= bases.size( ); j++ ) {
    base_t base2 = *b_iter; ++b_iter; // get next base
    // convert (base1, base2) into a color:
    colors.Set( j, SOLiD_toColor( base1, base2) );
    base1 = base2; // shift and prepare to retrieve the next base
  }
}

/// \brief Converts a collection of base space sequences into
/// a collection of colorspace representations of these sequences.
///
/// Colorspace representation of the sequence: first base of the 
/// sequence followed by colors, each representing a pair of adjacent
/// bases in the original sequence.
void ColorEncode( const vecbasevector& bases, vecbasevector& colors )
{   
  colors = bases;
  for (size_t i = 0; i < bases.size( ); i++ ) {    
    const basevector & B = bases[i];
    basevector & C = colors[i];
    ColorEncode( B , C );
  } 
}


/// \brief Converts a sequence from color space representation into
/// a base space representation. Colorspace representation is assumed
/// to have a base as its first element, followed by SOLiD colors,
/// representing adjacent base pairs. The base space representation will
/// start, accordingly, from the same base, followed by bases computed
/// by "unwinding" each preceding base and color back into a next base.
///
/// @param colors [in] input sequence in color space
/// @param bases [out] base space representation of the input sequence

void ColorDecode( const basevector & colors, basevector & bases ) {    
  bases.resize(colors.size());

  if ( colors.size() == 0 ) return ; // nothing to do 
  basevector::const_iterator c_iter = colors.Begin();
  base_t base = *c_iter; ++c_iter;
  
  bases.Set(0,base);

  for ( unsigned int j = 1 ; j < colors.size() ; j++ ) {
    SOLiD_color_t c = *c_iter; ++c_iter; // get next color
    // convert previous base and current color into next base:
    base = SOLiD_toBase(base,c); 
    bases.Set(j,base);
  }
}


/// \brief Converts a collection of color space sequences into
/// a collection of base space representations of these sequences.
///
/// Colorspace representation of the sequence must contain first base 
/// of the sequence followed by colors, each representing a pair of 
/// adjacent bases.

void ColorDecode( const vecbasevector& colors, vecbasevector& bases ) {    

  bases = colors;
  for (size_t i = 0; i < colors.size( ); i++ ) {    
    basevector & B = bases[i];
    const basevector & C = colors[i];
    ColorDecode( C , B );
  } 
}


/// Prints colorvector \c b into the output stream \c out adding VT220 
/// escape sequences to turn on colors. NOTE: as there is currently no difference
/// between internal representations of "true" basevectors and colorvectors,
/// this method \em assumes that the passed vector contains colors and prints
/// the data accordingly (as 0,1,2,3 rather than A,C,G,T). This method \em does
/// \em not expect the first element of the colorvector to be a base (followed by 
/// colors). If the colorvector used does contain a "full" colorspace sequence
/// (starting with a base), the first base should be removed manually (or 
/// argument \c start=1 should be used). This method does \em not print
/// newline character after the vector (think out << b; - no newline).

void PrintAsColorsWithEsc(ostream & out, const basevector & b, unsigned int start, int length) {
  if ( length == -1 ) length = b.size() - start ;
  basevector::const_iterator it =  b.Begin(start);
  for ( int i = 0 ; i < length ; i++ ) {
    base_t c = *it; ++it;
    if (i > 0 && i % 60 == 0 ) out << endl;	
    switch (c) {
    case 0:
      out << "\033[1;34m0\033[0m";
      break;
    case 1:
	//printf("1");      
      out << "\033[1;32m1\033[0m";
      break;
    case 2:
      //printf("2");      
      out << "\033[1;33m2\033[0m";
      break;
    case 3:
      //printf("3"); 
      out << "\033[1;31m3\033[0m";
      break;
    default:
      ForceAssert(0 == 1);
      //cout << (int)operator[](i) << endl;
    }
  }
}


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

void PrintAsColorsHTML(ostream & out, const basevector & b, unsigned int start, int length) {
  if ( length == -1 ) length = b.size() - start ;
  if ( length == 0 ) return;
  basevector::const_iterator it =  b.Begin(start);
  out << "<font face=\"Courier\">" << endl;
  base_t prev_c = 100;
  for ( int i = 0 ; i < length ; i++ ) {
    base_t c = *it; ++it;
    if (i > 0 && i % 60 == 0 ) out << endl;
    if ( c !=prev_c ) {
      if ( i!=0 ) out << "</font>";
      switch (c) {
      case 0:
	out << "<font color=\"Blue\">";
	break;
      case 1:
	//printf("1");      
	out << "<font color=\"LimeGreen\">";
	break;
      case 2:
	//printf("2");      
	out << "<font color=\"Orange\">";
	break;
      case 3:
	//printf("3"); 
	out << "<font color=\"Red\">";
	break;
      default:
	ForceAssert(0 == 1);
	//cout << (int)operator[](i) << endl;
      }
      prev_c=c;
    } // end if ( c!=prev_c ) - color is now switched
    out << c;
  }
  out << "</font>" << endl << "</font>" << endl;
}


/// Prints colorvector \c b into the output stream. 
/// NOTE: as there is currently no difference
/// between internal representations of "true" basevectors and colorvectors,
/// this method \em assumes that the passed vector contains colors and prints
/// the data accordingly (as 0,1,2,3 rather than A,C,G,T). This method \em does
/// \em not expect the first element of the colorvector to be a base (followed by 
/// colors). If the colorvector used does contain a "full" colorspace sequence
/// (starting with a base), the first base should be removed manually (or 
/// argument \c start=1 should be used). This method does \em not print
/// newline character after the vector (think out << b; - no newline).

void PrintAsColors(ostream & out, const basevector & b, unsigned int start, int length ) {
  if ( length == -1 ) length = b.size() - start ;
  basevector::const_iterator it = b.Begin(start);
  for ( int i = 0 ; i < length ; i++ ) {
    if ( i > 0 && i % 60 == 0 ) out << endl;
    out << (char)('0' + *it);
    ++it;
  }
}


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
/// newline character after the vector (think out << b; - no newline).
void PrintReverseAsColors(ostream & out, const basevector & b, unsigned int start, int length) {
  if ( length == -1 ) length = b.size() - start ;
  basevector::const_reverse_iterator it = b.RBegin(start);
  for ( int i = 0 ; i < length ; i++ ) {
    if ( i > 0 && i % 60 == 0 ) out << endl;
    out << (char)('0' + *it);
    ++it;
  }
}


/// Prints visual alignment of colorvector b. Alignment is specified by the
/// look_align argument and the basespace reference is passed in ref argument.
/// b is assumed to be in pure colorspace representation; if the first base is known
/// it should be passed separately in <first_base> argument  
void PrintVisualColorAlign(ostream & out, const basevector & b, 
			    const vecbasevector & ref,
			    const look_align & la, 
			   base_t first_base ) {

   basevector baseref; // stretch of the original ref the read aligns to
   baseref.SetToSubOf( ref[ la.target_id ], la.StartOnTarget(), la.query_length+1 );
   if ( la.Rc1() ) baseref.ReverseComplement();

   basevector colorref; // same stretch of the ref in color representation
   ColorEncode(baseref,colorref); // note: colorref keeps leading base!
   basevector::const_iterator t_it = colorref.Begin();
   basevector::const_iterator b_it = b.Begin();
   base_t crefBase = *t_it;
   ++t_it;
   if ( first_base != crefBase && first_base != 255 ) out << '*';
   else out << ' ';
   PrintMismatchString(out,b_it, t_it, la.query_length);
   out << endl;
   if ( first_base != 255 ) out << as_base(first_base) ;
   else out << ' ';
   PrintAsColors(out,b);
   out << endl;
   if ( first_base != 255 ) out << colorref.at(0) ;
   else out << ' ';
   PrintAsColors(out , colorref, 1, la.query_length );
   out << endl;
   baseref.PrintBases(out,0,la.query_length+1);
   out << endl;
}

/// Builds full colorspace representation (with the leading primer base \c first_base) of the target sequence
/// pointed to by the alignment \c la on the basespace reference \c ref (the sequence is also properly reverse 
/// complemented is the alignment is reverse). In other words, if you have a (full) colorspace read and its alignment,
/// this method will extract the reference sequence (in proper orientation) the read's DNA molecule purportedly
/// represents. The resulting colorspace sequence is
/// written into \c result (previous content, if any, is overwritten and the vector is resized)
void FullColorspaceTargetSeq(const look_align &la, const vecbasevector &ref, base_t first_base, basevector &result) {
       basevector bref ; // this is the actual reference sequence at the location,
                         // to which the read aligned
       bref.SetToSubOf( ref[la.TargetId()], la.StartOnTarget(), la.QueryLength()+1 );
       if ( la.IsQueryRC() ) {
	 bref.ReverseComplement(); // we keep the read as is and reverse complement reference seq
                                   // here if the alignment was in reverse direction...
       }
       // append primer base from the read, convert to full colorspace sequence:
       ColorEncode(first_base, bref, result); 
}

/// Applies simple one pass correction to the full colorspace reference \c read aligned as specified by \c la
/// to the reference \c ref [note: the reference must be provided in the original basespace representation]. The
/// error-corrected full colorspace sequence of the read will be written into \c result (previous content, if any,
/// is destroyed, and the vector is resized).
///
/// One-pass correction strategy is as follows: 1) every single mismatch flanked by matches is corrected;
/// 2) mismatch followed by a discordant mismatch is immediately corrected; 3) mismatch (1) followed by a concordant
/// mismatch (2) is assumed to be a SNP and never corrected, the position (2) is treated as match for the perposes
/// of correcting subsequent colors. Note that when this strategy is used, then, for instance, a sequence of
/// 3 mismatches (1)(2)(3) such that both (1),(2) and (2),(3) are concordant will be always called as a SNP at (1)(2)
/// (read colors (1)(2) will not be corrected) and a mismatch at (3) (i.e. color at (3) will get corrected).
///
/// This method returns a trio of counts: \c first is the number of corrected mismatches (those that did not end up in SNPs), 
/// \c second is the number of discordant pairs (\em adjacent pairs are counted, i.e. *** is 1 pair and **** is 2 pairs),
/// and \c third is the number of SNPs (concordant mismatches in two adjacent positions (pairs)) chosen in one pass as 
/// described above, so that total number of mismatching colors is \c first + 2*\c third.
trio<unsigned int, unsigned int,unsigned int> 
CorrectColorRead( const basevector & read, const vecbasevector & ref, 
		  const look_align & la, basevector & result ) {

       unsigned int n_mismatches = 0;
       unsigned int n_SNPs = 0; // number of concordant pairs
       unsigned int n_wrong_pairs  = 0; // number of discordant pairs
                                        // NOTE: these are not counted "combinatorially"
                                        // in a sense that '***' could be 2 (*overlapping*) pairs.
                                        // Instead, a continuous run of k mismatches that
                                        // do not form a single concordant pair (it would break the run!)
                                        // is counted as [k/2] *adjacent* discordant pairs

       result = read;

       unsigned int tig = la.TargetId();          // contig to which the current alignment was performed
       
       ForceAssertEq( la.QueryLength(), read.size() - 2 );

       // cref_conv_full sequence declared below will contain reference base sequence, 
       // prefixed with a primer base, and then converted into full colorspace representation.
       // Example: the full colorspace representation of the read is R=T012300 (where T is
       // the primer base), and its (meaningful) color sequence Rc=12300 aligns to (basespace) sequence
       // B=ACTAAA on the reference (check that it's 12300 indeed in pure colorspace representation!).
       // As this is where the read's actual DNA molecule purpotedly came from, we prepend the reference 
       // sequence B with read's primer base (T) and convert resulting
       // TACTAAA sequence into full colorspace representation: T312300. This way we can see that
       // mismatch at the first base posion (reference sequence starts with A, while read says T0=T)
       // is probably a sequencing error, not a SNP - because if it was a SNP, we would observe
       // 2 mismatches in the colorspace: if a true SNP were there and we'd have, say, TCTAAA on the ref, 
       // then prefixed with T this would give TTCTAAA -> T022300 in colorspace). 
       // Hence after prepending each aligned
       // sequence on the ref with a primer base we can treat the first base the same way as any other base:
       // only concordant color mismatch pairs constitute a SNP!

       basevector cref_conv_full;  
       FullColorspaceTargetSeq(la, ref, read[0], cref_conv_full);

       // we now have the color read and color-converted ref in full colorspace
       // representation. We can now match position by position (i.e. color by color)
       
       // we do not care about the firts base (primer base); set iterators to the first color:
       basevector::const_iterator read_iter = read.Begin(1); 
       basevector::const_iterator read_end = read.End();
       basevector::const_iterator cref_iter = cref_conv_full.Begin(1);
       
       // in the following, we scan in parallel the read and reference colorvectors 
       // keeping a state variable 'mismatch' that remembers if we have seen a (color) mismatch 
       // in the previous position.
       // The transitions between states (and correcting mismatches)
       // are governed by the following logic:
       // if we are in mismatch = false state ->
       //      if no mismatch in the current position -> stay in the same state
       //      if mismatch in the current position -> initiate delay: 
       //                                switch into mismatch=true state (i.e. wait until we see if
       //                                there is a mismatch in the following position and if it is
       //                                concordant so that we are actually seeing a SNP)
       // if we are in mismatch = true state ->
       //      if no mismatch in the current position -> mismatch in the previous position was
       //                              standalone (i.e. sequencing error); correct the read's color
       //                              in the previous position, switch to mismatch=false
       //      if mismatch in the current position -> check if mismatch in the current position agrees
       //                              with the mismatch in the previous position
       //             if the two mismatches agree -> it's a SNP.
       //                              switch to mismatch = false; DON'T correct the read.
       //             if mismatches do not agree -> correct the read's
       //                              color in the previous position, and stay
       //                              in mismatch = true state (i.e. we still have one-position delay and 
       //                              we are waiting -  we know only that the previous mismatch was wrong 
       //                              and we just corrected it, but we still
       //                              may receive a mismatch at the next step that will agree with the
       //                              current one and make it SNP-like!)

       // initializations:
       unsigned int pos = 1; // position on the corrected color read (if we are in 'mismatch' state, we may be
                             // 1 position behind where the read and ref iterators point to, see above!)
                             // we set 1 right away, because we do not care about the leading primer base!

       base_t read_color, ref_color; // will keep read and reference sequences' colors at the current position
       base_t read_color_prev = *read_iter; ++read_iter; // will keep the color we've seen at the previous position in the read
       base_t ref_color_prev = *cref_iter; ++cref_iter; // will keep the color we've seen at the previous position in the ref

       bool mismatch = false; 
       bool in_pair = false; // this is used to remember what happened *2* bases back, 
                             // so that we could correctly count wrong (discordant) pairs
       if ( read_color_prev != ref_color_prev ) mismatch = true; // mismatch at the very first color!
       else {
	 // if there's no mismatch then color at position pos in the corrected read is the same 
	 // as in the read. There's no delay, advance to the next position in the color read
	 pos++; 
       }

       // define z:=current position in the read and colorspace ref sequences (where the iterators point to)
       // precondition: z==1; (mismatch==false && pos==2 ) || (mismatch==true && pos==1)
       // precondition: read colors up to position pos-1 inclusive are verified/corrected

       while ( read_iter != read_end ) { // until the very end of the read/ref sequences

	   // the next two lines perfom z++, unconditionally. (z->z+1)
	   read_color = *read_iter; ++read_iter; // read next colors from read and ref
	   ref_color = *cref_iter; ++cref_iter;

	   // invariant: (it's more clear at this point in the loop):
	   //
	   //      read_color==read[z];        ref_color==cref_conv_full[z];
	   // read_color_prev==read[z-1]; ref_color_prev==cref_conv_full[z-1];
	   // ( mismatch==true && pos==z-1 ) || (mismatch==false && pos=z) 
	   // read colors up to position pos-1 are verified/corrected

	   if ( read_color == ref_color ) {
	     // no mismatch at current color position z
	     if ( mismatch ) { // did we have a mismatch at previous color position z-1?

	       // here mismatch == true => pos==z-1 (because of the invariant)

	       // correct color read at previous position, advance to current pos
	       result.Set(pos++, ref_color_prev);
	       mismatch = false;
	       n_mismatches++;
	       if ( in_pair ) {
		 // if the previous mismatch was the second one in a pair, count wrong pair:
		 in_pair=false;
		 n_wrong_pairs++;
	       }
	       // now ( mismatch==false && pos==z) ; read colors up to pos verified/corrected

	     }
	     else { // EMPTY ON PURPOSE - nothing to do

	       // here mismatch==false => pos==z; unchanged in this else-clause
	       // read colors up to pos are verified/corrected

	     }
	     // no mismatch in the current position, advance to the next:
	     pos++;
	     // now: mismatch==false; pos==z+1; read colors up to pos-1 are verified/corrected
	   } else {
	     // oops, we got mismatch:
	     if ( mismatch ) { 

	       // if we already had a mismatch at previous position, then we got 2 mismatches in a row!
	       
	       // here mismatch==true => pos==z-1

	       if ( (read_color_prev ^ read_color) == (ref_color_prev ^ ref_color) ) { // do we have a SNP?
		 mismatch = false;
		 in_pair=false; // even if the previous mismatch was the second one in a pair,
		                // we count this mismatch into a SNP here, so there's no preceding pair!
		                // reset flag and do not count wrong pair
		 n_SNPs++;
		 // on SNPs we don't correct the colors in the edited read; just advance
		 // past these two positions:
		 pos++; 
		 pos++;
		 
		 // now mismatch==false; pos==z+1
		 // read colors up to pos-1 are verified/corrected

	       } else { // we got second mismatch in a row, but it's discordant - not a SNP!
		 n_mismatches++;
		 if ( in_pair ) {
		   // we already had two mismatches in a row prior to this mismatch
		   in_pair=false;
		   n_wrong_pairs++;
		 } else {
		   // this one is a second mismatch in row, but no mismatches 2 colors back
		   in_pair = true ; // mark, but not count yet: mismatch we are looking at may still be first one in a SNP
		 }
		 // previous mismatch was incorrect (misread), correct the color read in that
		 // position and advance to the current. Leave the current position for later decision.
		 result.Set(pos++,ref_color_prev);
		 // now mismatch == true; pos==z;
		 // read colors up to pos-1 are verified/corrected

	       }
	     } else {

	       // here mismatch==false => pos==z
	       mismatch = true; // we just got a mismatch after a match... change state, wait and see...
	       // now mismatch==true; pos==z
	       // read colors up to pos-1 are verified/corrected

	     }
	   }
	   read_color_prev = read_color; 
	   ref_color_prev = ref_color;

	   // now read_color_prev==read[z]; ref_color_prev==cref_conv_full[z];
	   // (mismatch==true && pos==z) || ( mismatch==false && pos==z+1)
	   // read colors up to pos-1 inclusive are verified/corrected
	   // after the first two lines at the beginning of the loop, the invariant will hold!

       }  // while ( read_iter!=read_end ) - end of the loop over all positions in color read
	 
       // post-processing

       if ( mismatch ) { 
	 // here mismatch==true => pos=z-1; z==read.length - 1
	 // if we ended up waiting with an unpaired mismatch at last position 
	 // (all other cases are taken care of!)
	 n_mismatches++;
	 if ( in_pair ) n_wrong_pairs++;
	 result.Set(pos++,ref_color_prev); // and correct it
       }

       return(make_trio(n_mismatches, n_wrong_pairs, n_SNPs));
}

