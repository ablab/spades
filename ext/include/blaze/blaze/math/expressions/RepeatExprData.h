//=================================================================================================
/*!
//  \file blaze/math/expressions/RepeatExprData.h
//  \brief Header file for the implementation of the RepeatExprData class template
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

#ifndef _BLAZE_MATH_EXPRESSIONS_REPEATEXPRDATA_H_
#define _BLAZE_MATH_EXPRESSIONS_REPEATEXPRDATA_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/StaticAssert.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  FORWARD DECLARATIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class template for the data members of repeater expression classes.
// \ingroup math
//
// The auxiliary RepeatExprData class template represents an abstraction of the data members of
// repeater expression template classes. The necessary set of data member is selected depending
// on the number of compile time repeater arguments.
*/
template< size_t Dim        // Number of dimensions
        , size_t... CRAs >  // Compile time repeater arguments
class RepeatExprData;
//*************************************************************************************************




//=================================================================================================
//
//  BASE TEMPLATE FOR MANY COMPILE TIME REPEATER ARGUMENTS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Base template of the RepeatExprData class template.
// \ingroup math
//
// The base template of RepeatExprData adapts the class template to the requirements of a many
// compile time repeater arguments.
*/
template< size_t Dim      // Number of dimensions
        , size_t... Rs >  // Compile time repetitions
class RepeatExprData
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   constexpr RepeatExprData() noexcept;

   RepeatExprData( const RepeatExprData& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~RepeatExprData() = default;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   RepeatExprData& operator=( const RepeatExprData& ) = delete;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   template< size_t I >
   static constexpr size_t repetitions() noexcept;
   //@}
   //**********************************************************************************************

 private:
   //**Compile time checks*************************************************************************
   BLAZE_STATIC_ASSERT( Dim == sizeof...( Rs ) );
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief The constructor for RepeatExprData.
*/
template< size_t Dim      // Number of dimensions
        , size_t... Rs >  // Compile time repetitions
constexpr RepeatExprData<Dim,Rs...>::RepeatExprData() noexcept
{}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of repetitions of the repeater expression for the given dimension.
//
// \return The number of repetitions for the given dimension.
*/
template< size_t Dim      // Number of dimensions
        , size_t... Rs >  // Compile time repetitions
template< size_t I >      // Dimension index
constexpr size_t RepeatExprData<Dim,Rs...>::repetitions() noexcept
{
   BLAZE_STATIC_ASSERT( I < Dim );
   constexpr size_t reps[] = { Rs... };
   return reps[I];
}

/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR ZERO COMPILE TIME REPEATER ARGUMENTS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the RepeatExprData class template for zero compile time repeater
//        arguments.
// \ingroup math
//
// This specialization of RepeatExprData adapts the class template to the requirements of zero
// compile time repeater arguments.
*/
template< size_t Dim >  // Number of dimensions
class RepeatExprData<Dim>
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   template< typename... Reps >
   inline RepeatExprData( Reps... reps ) noexcept;

   RepeatExprData( const RepeatExprData& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~RepeatExprData() = default;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   RepeatExprData& operator=( const RepeatExprData& ) = delete;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   template< size_t I >
   inline size_t repetitions() const noexcept;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   const size_t repetitions_[Dim];  //!< The number of repetitions of the repeater expression.
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief The constructor for RepeatExprData.
//
// \param reps The number of repetitions of the repeater expression.
*/
template< size_t Dim >        // Number of dimensions
template< typename... Reps >  // Runtime repetitions
inline RepeatExprData<Dim>::RepeatExprData( Reps... reps ) noexcept
   : repetitions_{ reps... }  // The number of repetitions of the repeater expression
{
   BLAZE_STATIC_ASSERT( sizeof...( Reps ) == Dim );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of repetitions of the repeater expression for the given dimension.
//
// \return The number of repetitions for the given dimension.
*/
template< size_t Dim >  // Number of dimensions
template< size_t I >    // Dimension index
inline size_t RepeatExprData<Dim>::repetitions() const noexcept
{
   BLAZE_STATIC_ASSERT( I < 3UL );
   return repetitions_[I];
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
