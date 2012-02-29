///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file Exit.h
 * \author tsharpe
 * \date Feb 6, 2009
 *
 * \brief A replacement for the ::exit(int) function.
 *
 * Comes in two flavors:
 * In a single-threaded program it calls ::exit(0) if the arg is zero,
 * and calls abort() otherwise.
 * In a multi-threaded program, it calls ::exit(0) if the arg is zero,
 * and call abort() otherwise for the main thread only.  For all other
 * threads it call pthread_exit.
 */
#ifndef SYSTEM_EXIT_H_
#define SYSTEM_EXIT_H_

namespace CRD
{
void exit( int ) __attribute__((__noreturn__));
} // end namespace CRD

inline void TracebackThisProcess() { CRD::exit(1); }

#endif /* SYSTEM_EXIT_H_ */
