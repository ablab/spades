//=================================================================================================
/*!
//  \file blaze/config/Debugging.h
//  \brief Configuration of the debugging policy of the Blaze library
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


//*************************************************************************************************
/*!\brief Compilation switch for the (de-)activation of the debug mode.
// \ingroup config
//
// This compilation switch enables/disables the debug mode of the Blaze library. In case the
// switch is set to 1 (i.e. in case the debug mode is enabled), the Blaze library is allowed
// to perform additional runtime checks and to modify user-defined settings to guarantee a
// full and thorough debugging process. In case the switch is set to 0 (i.e. the debug mode
// is disabled), no additional checks are added and no settings are modified.
//
// Possible settings for the debug mode switch:
//  - Deactivated: \b 0 (default)
//  - Activated  : \b 1
//
// \note It is possible to (de-)activate the debug mode via command line or by defining this
// symbol manually before including any Blaze header file:

   \code
   g++ ... -DBLAZE_USE_DEBUG_MODE=1 ...
   \endcode

   \code
   #define BLAZE_USE_DEBUG_MODE 1
   #include <blaze/Blaze.h>
   \endcode
*/
#ifndef BLAZE_USE_DEBUG_MODE
#define BLAZE_USE_DEBUG_MODE 0
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for function traces.
// \ingroup config
//
// This compilation switch triggers the use of function traces. In case the switch is set to
// 1, function traces via the BLAZE_FUNCTION_TRACE are enabled and trace information is written
// to the console via \c std::cerr.
//
// Possible settings for the function trace switch:
//  - Deactivated: \b 0 (default)
//  - Activated  : \b 1
//
// \note It is possible to (de-)activate function traces via command line or by defining this
// symbol manually before including any Blaze header file:

   \code
   g++ ... -DBLAZE_USE_FUNCTION_TRACES=1 ...
   \endcode

   \code
   #define BLAZE_USE_FUNCTION_TRACES 1
   #include <blaze/Blaze.h>
   \endcode
*/
#ifndef BLAZE_USE_FUNCTION_TRACES
#define BLAZE_USE_FUNCTION_TRACES 0
#endif
//*************************************************************************************************
