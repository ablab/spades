//=================================================================================================
/*!
//  \file blaze/math/views/column/ColumnData.h
//  \brief Header file for the implementation of the ColumnData class template
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

#ifndef _BLAZE_MATH_VIEWS_COLUMN_COLUMNDATA_H_
#define _BLAZE_MATH_VIEWS_COLUMN_COLUMNDATA_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/MaybeUnused.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class template for the data members of the Column class.
// \ingroup column
//
// The auxiliary ColumnData class template represents an abstraction of the data members of
// the Column class template. The necessary set of data members is selected depending on the
// number of compile time column arguments.
*/
template< size_t... CCAs >  // Compile time column arguments
class ColumnData
{};
//*************************************************************************************************




//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR ZERO COMPILE TIME COLUMN INDICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the ColumnData class template for zero compile time column arguments.
// \ingroup column
//
// This specialization of ColumnData adapts the class template to the requirements of zero compile
// time column arguments.
*/
template<>
class ColumnData<>
{
 public:
   //**Compile time flags**************************************************************************
   //! Compilation flag for compile time optimization.
   /*! The \a compileTimeArgs compilation flag indicates whether the view has been created by
       means of compile time arguments and whether these arguments can be queried at compile
       time. In that case, the \a compileTimeArgs compilation flag is set to \a true, otherwise
       it is set to \a false. */
   static constexpr bool compileTimeArgs = false;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   template< typename... RCAs >
   explicit inline ColumnData( size_t index, RCAs... args );

   ColumnData( const ColumnData& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~ColumnData() = default;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   ColumnData& operator=( const ColumnData& ) = delete;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline size_t column() const noexcept;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   const size_t column_;  //!< The index of the column in the matrix.
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief The constructor for ColumnData.
//
// \param index The index of the column.
// \param args The optional column arguments.
*/
template< typename... RCAs >  // Optional column arguments
inline ColumnData<>::ColumnData( size_t index, RCAs... args )
   : column_( index )  // The index of the column in the matrix
{
   MAYBE_UNUSED( args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the index of the column of the underlying dense matrix.
//
// \return The index of the column.
*/
inline size_t ColumnData<>::column() const noexcept
{
   return column_;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR ONE COMPILE TIME COLUMN INDEX
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the ColumnData class template for a single compile time column argument.
// \ingroup column
//
// This specialization of ColumnData adapts the class template to the requirements of a single
// compile time column argument.
*/
template< size_t I >   // Compile time column index
class ColumnData<I>
{
 public:
   //**Compile time flags**************************************************************************
   //! Compilation flag for compile time optimization.
   /*! The \a compileTimeArgs compilation flag indicates whether the view has been created by
       means of compile time arguments and whether these arguments can be queried at compile
       time. In that case, the \a compileTimeArgs compilation flag is set to \a true, otherwise
       it is set to \a false. */
   static constexpr bool compileTimeArgs = true;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   template< typename... RCAs >
   explicit inline ColumnData( RCAs... args );

   ColumnData( const ColumnData& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~ColumnData() = default;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   ColumnData& operator=( const ColumnData& ) = delete;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   static constexpr size_t column() noexcept;
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief The constructor for ColumnData.
//
// \param args The optional column arguments.
*/
template< size_t I >          // Compile time column index
template< typename... RCAs >  // Optional column arguments
inline ColumnData<I>::ColumnData( RCAs... args )
{
   MAYBE_UNUSED( args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the index of the column of the underlying dense matrix.
//
// \return The index of the column.
*/
template< size_t I >  // Compile time column index
constexpr size_t ColumnData<I>::column() noexcept
{
   return I;
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
