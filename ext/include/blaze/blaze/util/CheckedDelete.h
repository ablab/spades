//=================================================================================================
/*!
//  \file blaze/util/CheckedDelete.h
//  \brief Type-checked delete operations
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

#ifndef _BLAZE_UTIL_CHECKEDDELETE_H_
#define _BLAZE_UTIL_CHECKEDDELETE_H_


namespace blaze {

//=================================================================================================
//
//  CHECKED DELETE OPERATIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Checked delete operations */
//@{
template< typename T > void checkedDelete( T* ptr );
template< typename T > void checkedArrayDelete( T* ptr );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Type-checked \c delete operation.
// \ingroup util
//
// \param ptr The pointer to be deleted.
// \return void
//
// This function frees the given pointer resource via \c delete. Note that the \c delete operation
// is NOT permitted for incomplete types (i.e. declared but undefined data types). The attempt
// to use this function for a pointer to an object of incomplete type results in a compile time
// error!
*/
template< typename T >
void checkedDelete( T* ptr )
{
   using TypeMustBeComplete = char[ sizeof(T)? 1 : -1 ];
   (void) sizeof(TypeMustBeComplete);
   delete ptr;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Type-checked \c delete[] operation.
// \ingroup util
//
// \param ptr The pointer to the array to be deleted.
// \return void
//
// This function frees the given pointer resource via \c delete[]. Note that the \c delete[]
// operation is NOT permitted for incomplete types (i.e. declared but undefined data types).
// The attempt to use this function for a pointer to an array of objects of incomplete type
// results in a compile time error!
*/
template< typename T >
void checkedArrayDelete( T* ptr )
{
   using TypeMustBeComplete = char[ sizeof(T)? 1 : -1 ];
   (void) sizeof(TypeMustBeComplete);
   delete[] ptr;
}
//*************************************************************************************************

} // namespace blaze

#endif
