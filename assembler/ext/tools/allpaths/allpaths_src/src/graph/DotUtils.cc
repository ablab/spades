/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "graph/DotUtils.h"
#include "random/Random.h"

vec< vec<dot_color_t> > DotUtils::dotColors_;

static struct DotUtilsInitializer {
  DotUtilsInitializer() {
    DotUtils::InitDotColors();
  }
} myIniter;

void DotUtils::InitDotColors() {
  const char *whites[] = {
    "antiquewhite", "azure", "bisque",
    "blanchedalmond",
    "cornsilk", "floralwhite", "gainsboro",
    "ghostwhite", "honeydew", "ivory", "lavender",
    "lavenderblush", "lemonchiffon",
    "linen", "mintcream", "mistyrose",
    "moccasin", "navajowhite", "oldlace",
    "papayawhip", "peachpuff", "seashell",
    "snow", "thistle", "wheat", "white", "whitesmoke"
    };
  const char *greys[] = {
    "darkslategray", "dimgray", "gray", "lightgray",
    "lightslategray", "slategray"
  };
  const char *blacks[] = {
    "black"
  };
  const char *reds[] = {
    "coral", "crimson", "darksalmon",
    "deeppink", "firebrick", "hotpink",
    "indianred", "lightpink", "lightsalmon",
    "maroon", "mediumvioletred",
    "orangered", "palevioletred",
    "pink", "red", "salmon", "tomato", "violetred"
  };
  const char *browns[] = {
    "beige", "brown", "burlywood", "chocolate",
    "darkkhaki", "khaki", "peru", "rosybrown",
    "saddlebrown", "sandybrown", "sienna", "tan"
  };
  const char *oranges[] = {
    "darkorange", "orange", "orangered"
  };
  const char *yellows[] = {
    "darkgoldenrod", "gold", "goldenrod",
    "greenyellow", "lightgoldenrod", "lightgoldenrodyellow",
    "lightyellow", "palegoldenrod", "yellow", "yellowgreen"
  };
  const char *greens[] = {
    "chartreuse", "darkgreen", "darkolivegreen", "darkseagreen",
    "forestgreen", "green", "greenyellow", "lawngreen", "lightseagreen",
    "limegreen", "mediumseagreen", "mediumspringgreen",
    "mintcream", "olivedrab", "palegreen", "seagreen", "springgreen",
    "yellowgreen"
  };
  const char *cyans[] = {
    "aquamarine", "cyan", "darkturquoise", "lightcyan", "mediumaquamarine",
    "mediumturquoise", "paleturquoise", "turquoise"
  };
  const char *blues[] = {
    "aliceblue", "blue", "blueviolet", "cadetblue",
    "cornflowerblue", "darkslateblue", "deepskyblue", "dodgerblue",
    "indigo", "lightblue", "lightskyblue", "lightslateblue",
    "mediumblue", "mediumslateblue", "midnightblue",
    "navy", "navyblue", "powderblue", "royalblue", "skyblue",
    "slateblue", "steelblue"
  };
  const char *magentas[] = {
    "blueviolet", "darkorchid", "darkviolet",
    "magenta", "mediumorchid", "mediumpurple", "mediumvioletred",
    "orchid", "palevioletred", "plum", "purple", "violet", "violetred"
  };
  //AddDotColors( whites, sizeof(whites) );
  AddDotColors( greys, sizeof(greys) );
  AddDotColors( blacks, sizeof(blacks) );
  AddDotColors( reds, sizeof(reds) );
  AddDotColors( oranges, sizeof(oranges) );
  AddDotColors( yellows, sizeof(yellows) );
  AddDotColors( greens, sizeof(greens) );
  AddDotColors( cyans, sizeof(cyans) );
  AddDotColors( blues, sizeof(blues) );
  AddDotColors( magentas, sizeof(magentas) );
}

void DotUtils::AddDotColors( const char *colors[], int ncolors ) {
  dotColors_.resize(dotColors_.size()+1);
  for ( unsigned int i = 0; i < ( ncolors / sizeof(const char *) ); i++ )
    dotColors_.back().push_back( dot_color_t( colors[i] ) );
}

dot_color_t DotUtils::RandomColor() {
#if 0  
  int colorGroup = randomx() % dotColors_.size();
  int colorId = randomx() % ( dotColors_[colorGroup].size() ) ;
  return dotColors_[ colorGroup  ][ colorId ];
#endif
  char clr[1024];
  sprintf(clr, "#%02x%02x%02x", unsigned(randomx() % 0xff), unsigned(randomx() % 0xff),
	  unsigned(randomx() % 0xff) );
  return dot_color_t( clr );
}

void DotUtils::RandomFgBg( dot_color_t& fg, dot_color_t& bg ) {
  unsigned r = randomx() % 0xff, g = randomx() % 0xff, b = randomx() % 0xff;
  char clr[1024];
  sprintf( clr, "#%02x%02x%02x", r, g, b );
  fg = clr;
  sprintf( clr, "#%02x%02x%02x", 255-r, 255-g, 255-b );
  bg = clr;
}




