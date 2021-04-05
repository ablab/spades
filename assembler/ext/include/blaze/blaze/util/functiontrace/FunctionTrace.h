//=================================================================================================
/*!
//  \file blaze/util/functiontrace/FunctionTrace.h
//  \brief Header file for the FunctionTrace class
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

#ifndef _BLAZE_UTIL_TRACING_FUNCTIONTRACE_H_
#define _BLAZE_UTIL_TRACING_FUNCTIONTRACE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#if BLAZE_HPX_PARALLEL_MODE
#  include <hpx/include/threads.hpp>
#elif BLAZE_CPP_THREADS_PARALLEL_MODE
#  include <thread>
#elif BLAZE_BOOST_THREADS_PARALLEL_MODE
#  include <boost/thread/thread.hpp>
#elif BLAZE_OPENMP_PARALLEL_MODE
#  include <omp.h>
#endif

#include <iostream>
#include <new>
#include <sstream>
#include <string>
#include <blaze/system/Debugging.h>
#include <blaze/system/Signature.h>
#include <blaze/system/SMP.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief RAII object for function tracing.
// \ingroup util
//
// The FunctionTrace class is an auxiliary helper class for the tracing of function calls. It
// is implemented as a wrapper around \c std::cerr and is responsible for the atomicity of the
// output of trace information.
*/
class FunctionTrace
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   inline FunctionTrace( const std::string& file, const std::string& function );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   inline ~FunctionTrace();
   //@}
   //**********************************************************************************************

   //**Forbidden operations************************************************************************
   /*!\name Forbidden operations */
   //@{
   FunctionTrace( const FunctionTrace& ) = delete;
   FunctionTrace( FunctionTrace&& ) = delete;

   FunctionTrace& operator=( const FunctionTrace& ) = delete;
   FunctionTrace& operator=( FunctionTrace&& ) = delete;

   void* operator new  ( std::size_t ) = delete;
   void* operator new[]( std::size_t ) = delete;
   void* operator new  ( std::size_t, const std::nothrow_t& ) noexcept = delete;
   void* operator new[]( std::size_t, const std::nothrow_t& ) noexcept = delete;

   void operator delete  ( void* ) noexcept = delete;
   void operator delete[]( void* ) noexcept = delete;
   void operator delete  ( void*, const std::nothrow_t& ) noexcept = delete;
   void operator delete[]( void*, const std::nothrow_t& ) noexcept = delete;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   std::string file_;      //!< The file name the traced function is contained in.
   std::string function_;  //!< The name of the traced function.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the FunctionTrace class.
//
// \param file The name of the file the traced function is contained in
// \param function The name of the traced function
*/
inline FunctionTrace::FunctionTrace( const std::string& file, const std::string& function )
   : file_    ( file     )  // The file name the traced function is contained in
   , function_( function )  // The name of the traced function
{
   std::ostringstream oss;
   oss << " + ";

#if BLAZE_HPX_PARALLEL_MODE
   oss << "[Thread " << hpx::this_thread::get_id() << "]";
#elif BLAZE_CPP_THREADS_PARALLEL_MODE
   oss << "[Thread " << std::this_thread::get_id() << "]";
#elif BLAZE_BOOST_THREADS_PARALLEL_MODE
   oss << "[Thread " << boost::this_thread::get_id() << "]";
#elif BLAZE_OPENMP_PARALLEL_MODE
   oss << "[Thread " << omp_get_thread_num() << "]";
#endif

   oss << " Entering function '" << function_ << "' in file '" << file_ << "'\n";
   std::cerr << oss.str();
}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Destructor for the FunctionTrace class.
*/
inline FunctionTrace::~FunctionTrace()
{
   std::ostringstream oss;
   oss << " - ";

#if BLAZE_OPENMP_PARALLEL_MODE
   oss << "[Thread " << omp_get_thread_num() << "]";
#elif BLAZE_CPP_THREADS_PARALLEL_MODE
   oss << "[Thread " << std::this_thread::get_id() << "]";
#elif BLAZE_BOOST_THREADS_PARALLEL_MODE
   oss << "[Thread " << boost::this_thread::get_id() << "]";
#elif BLAZE_HPX_PARALLEL_MODE
   oss << "[Thread " << hpx::this_thread::get_id() << "]";
#endif

   oss << " Leaving function '" << function_ << "' in file '" << file_ << "'\n";
   std::cerr << oss.str();
}
//*************************************************************************************************

} // namespace blaze

#endif
