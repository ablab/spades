//=================================================================================================
/*!
//  \file blaze/math/expressions/ExpandExprData.h
//  \brief Header file for the implementation of the ExpandExprData class template
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

#ifndef _BLAZE_MATH_EXPRESSIONS_EXPANDEXPRDATA_H_
#define _BLAZE_MATH_EXPRESSIONS_EXPANDEXPRDATA_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class template for the data members of expansion expression classes.
// \ingroup math
//
// The auxiliary ExpandExprData class template represents an abstraction of the data members of
// expansion expression template classes. The necessary set of data member is selected depending
// on the number of compile time expansion arguments.
*/
template< size_t... CEAs >  // Compile time expansion arguments
class ExpandExprData
{};
//*************************************************************************************************




//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR ZERO COMPILE TIME EXPANSION ARGUMENTS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the ExpandExprData class template for zero compile time expansion
//        arguments.
// \ingroup math
//
// This specialization of ExpandExprData adapts the class template to the requirements of zero
// compile time expansion arguments.
*/
template<>
class ExpandExprData<>
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   inline ExpandExprData( size_t expansion ) noexcept;

   ExpandExprData( const ExpandExprData& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~ExpandExprData() = default;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   ExpandExprData& operator=( const ExpandExprData& ) = delete;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline size_t expansion() const noexcept;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   const size_t expansion_;  //!< The expansion of the expansion expression.
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief The constructor for ExpandExprData.
//
// \param expansion The expansion of the expansion expression.
*/
inline ExpandExprData<>::ExpandExprData( size_t expansion ) noexcept
   : expansion_( expansion )  // The expansion of the expansion expression
{}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the expansion of the expansion expression.
//
// \return The expansion.
*/
inline size_t ExpandExprData<>::expansion() const noexcept
{
   return expansion_;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR ONE COMPILE TIME EXPANSION ARGUMENT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the ExpandExprData class template for one compile time expansion
//        argument.
// \ingroup math
//
// This specialization of ExpandExprData adapts the class template to the requirements of a
// single compile time expansion argument.
*/
template< size_t E >  // Compile time expansion
class ExpandExprData<E>
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   constexpr ExpandExprData() noexcept;

   ExpandExprData( const ExpandExprData& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~ExpandExprData() = default;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   ExpandExprData& operator=( const ExpandExprData& ) = delete;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   static constexpr size_t expansion() noexcept;
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief The constructor for ExpandExprData.
*/
template< size_t E >  // Compile time expansion
constexpr ExpandExprData<E>::ExpandExprData() noexcept
{}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the expansion of the expansion expression.
//
// \return The expansion.
*/
template< size_t E >  // Compile time expansion
constexpr size_t ExpandExprData<E>::expansion() noexcept
{
   return E;
}

/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
