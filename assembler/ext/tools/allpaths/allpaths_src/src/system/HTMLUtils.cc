/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "system/HTMLUtils.h"

String HTMLHeadFrames( String title, String scripts, String bodyProps ) {
 return
   //  "<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01 Frameset//EN\" \"http://www.w3.org/TR/html4/frameset.dtd\">"
   "<html><head>"
   "<meta http-equiv=\"Content-Type\" content=\"text/html; charset=UTF-8\">"
   "<title>" + title + "</title>"
   "<META HTTP-EQUIV=\"Expires\" CONTENT=\"Thu, 1 June 2000 23:59:00 GMT\">"
   "<META HTTP-EQUIV=\"pragma\" CONTENT=\"no-cache\">" + scripts + "</head>\n\n";
}

String HTMLTailFrames() {
  return "\n</html>\n";
}

 
String HTMLHead( String title, String scripts, String bodyProps ) {
 return
   //   "<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\"  \"http://www.w3.org/TR/html4/loose.dtd\">"
   "<html><head>"
   "<meta http-equiv=\"Content-Type\" content=\"text/html; charset=UTF-8\">"
   "<title>" + title + "</title>"
   "<META HTTP-EQUIV=\"Expires\" CONTENT=\"Thu, 1 June 2000 23:59:00 GMT\">"
   "<META HTTP-EQUIV=\"pragma\" CONTENT=\"no-cache\">" + scripts + "</head>\n<body "
   + bodyProps +
   ">\n\n";
}

String HTMLTail() {
  return "\n</body></html>\n";
}

 

