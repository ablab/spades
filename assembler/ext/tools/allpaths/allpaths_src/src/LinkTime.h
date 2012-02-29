///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file LinkTime.h
 * \author tsharpe
 * \date Oct 13, 2009
 *
 * \brief
 */
#ifndef LINKTIME_H_
#define LINKTIME_H_

// These declarations are defined in module LinkTimestamp.cc
//
// We've given this file a slightly different name (LinkTime.h vs.
// LinkTimestamp.cc) on purpose:
// We always recompile LinkTimestamp.cc at the start of the make process, and
// so we don't want any program to depend on it explicitly or it would force
// that program to be re-linked every time.

extern char const* OS_RELEASE;
extern char const* ARACHNE_RELEASE;
extern char const* LINK_TIMESTAMP;
extern char const* SVN_REVISION;

#endif /* LINKTIME_H_ */
