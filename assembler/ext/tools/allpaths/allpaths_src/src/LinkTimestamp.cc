///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#define STRINGIZE(X) #X
#define STRINGIZER(X) STRINGIZE(X)

char const* OS_RELEASE = STRINGIZER(MAKE_OS_RELEASE);
char const* ARACHNE_RELEASE = STRINGIZER(MAKE_RELEASE);
char const* LINK_TIMESTAMP = __DATE__ " " __TIME__;
char const* SVN_REVISION = STRINGIZER(SVN_VERSION);

char const* ThisCanBeReadUsingTheStringsProgram =
 "OS Release=" STRINGIZER(MAKE_OS_RELEASE) \
 ", Release=" STRINGIZER(MAKE_RELEASE) \
 ", Link Timestamp=" __DATE__ " " __TIME__ \
 ", SVN Revision=" STRINGIZER(SVN_VERSION);
