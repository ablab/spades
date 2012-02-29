/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "Vec.h"
#include "String.h"

/**
   Header: DotUtils.h

   Utils for working with .dot graphs.
*/

// Type: dot_color_t
// A string representing a dot color.
typedef String dot_color_t;

/**
   Class: DotUtils

   Various static methods to help with generating dot graphs.
*/
class DotUtils {
 public:
  // Method: DotColors
  // Returns a constant array of dot color names, grouped into groups where colors
  // in the same group are similar while colors in different groups are different.
  static const vec< vec<dot_color_t> >& DotColors() { if ( dotColors_.empty() ) InitDotColors(); return dotColors_; }

  static dot_color_t RandomColor();

  static void RandomFgBg( dot_color_t& fg, dot_color_t& bg );

 private:
  static vec< vec<dot_color_t> > dotColors_;

  static void InitDotColors();
  static void AddDotColors( const char *colors[], int ncolors );

  friend struct DotUtilsInitializer;
  
};

