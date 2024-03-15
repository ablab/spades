//=================================================================================================
/*!
//  \file blaze/util/FunctionTrace.h
//  \brief Header file for the function trace functionality
//
//  Copyright (C) 2012-2020 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. You can redistribute it and/or modify it under
//  the terms of the New (Revised) BSD License. Redistribution and use in source and binary
//  forms, with or without modification, are permitted provided that the following conditions
//  are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice, this list
//     of conditions and the following disclaimer in the documentation and/or other materials
//     provided with the distribution.
//  3. Neither the names of the Blaze development group nor the names of its contributors
//     may be used to endorse or promote products derived from this software without specific
//     prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
//  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
//  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
//  SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
//  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
//  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
//  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
//  DAMAGE.
*/
//=================================================================================================

#ifndef _BLAZE_UTIL_FUNCTIONTRACE_H_
#define _BLAZE_UTIL_FUNCTIONTRACE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/system/Debugging.h>
#if BLAZE_USE_FUNCTION_TRACES
#  include <blaze/util/functiontrace/FunctionTrace.h>
#endif




//=================================================================================================
//
//  BLAZE_FUNCTION_TRACE MACRO
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Function trace macro.
// \ingroup util
//
// This macro can be used to reliably trace function calls. In case function tracing is
// activated, function traces are written to the console via \c std::cerr. The following
// short example demonstrates how the function trace macro is used:

   \code
   int main( int argc, char** argv )
   {
      BLAZE_FUNCTION_TRACE;

      // ...
   }
   \endcode

// The macro should be used as the very first statement inside the function in order to
// guarantee that writing the function trace is the very first and last action of the
// function call.\n
// Function tracing can be enabled or disabled via the BLAZE_USE_FUNCTION_TRACES macro.
// If function tracing is activated, trace information of the following form will be written
// to \c std::cerr:

   \code
   + [Thread 0] Entering function 'int main()' in file 'TraceDemo.cpp'
   - [Thread 0] Leaving function 'int main()' in file 'TraceDemo.cpp'
   \endcode

// In case function tracing is deactivated, all function trace functionality is completely
// removed from the code, i.e. no function traces are logged and no overhead results from
// the BLAZE_FUNCTION_TRACE macro.
*/
#if BLAZE_USE_FUNCTION_TRACES
#  define BLAZE_FUNCTION_TRACE \
   blaze::FunctionTrace BLAZE_FUNCTION_TRACE_OBJECT( __FILE__, BLAZE_SIGNATURE )
#else
#  define BLAZE_FUNCTION_TRACE
#endif
//*************************************************************************************************

#endif
