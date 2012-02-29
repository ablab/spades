/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/**
   Header: HTMLUtils.h

   Utilities that help generate HTML pages.  Mostly used for presenting debugging
   output in html.
*/

#ifndef __INCLUDE_system_HTMLUtils_h
#define __INCLUDE_system_HTMLUtils_h

#include <iostream>
#include "String.h"
#include "SemanticTypes.h"

// Semantic Type: html_t
// A string of well-formed HTML.   This need not be a standalone
// HTML element, but should be a well-formed html in that all
// the elements are balanced and closed.
SemanticType( String, html_t );

// Semantic Type: html_attrs_t
// A string of well-formed HTML attributes, of them form A=1 B=2....  
SemanticType( String, html_attrs_t );

// Semantic Type: url_t
// A string representing a URL
typedef String url_t;

html_t HTMLHead( String title, html_t scripts = "", String bodyProps = "" );
html_t HTMLTail();
html_t HTMLHeadFrames( String title, html_t scripts = "", String bodyProps="" );
html_t HTMLTailFrames();

#endif
// #ifndef __INCLUDE_system_HTMLUtils_h
