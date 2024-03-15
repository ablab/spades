//=================================================================================================
/*!
//  \file blaze/math/views/band/BandData.h
//  \brief Header file for the implementation of the BandData class template
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

#ifndef _BLAZE_MATH_VIEWS_BAND_BANDDATA_H_
#define _BLAZE_MATH_VIEWS_BAND_BANDDATA_H_


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
/*!\brief Auxiliary class template for the data members of the Band class.
// \ingroup band
//
// The auxiliary BandData class template represents an abstraction of the data members of the
// Band class template. The necessary set of data member is selected depending on the number
// of compile time band arguments.
*/
template< ptrdiff_t... CBAs >  // Compile time band arguments
class BandData
{};
//*************************************************************************************************




//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR ZERO COMPILE TIME BAND INDICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the BandData class template for zero compile time band arguments.
// \ingroup band
//
// This specialization of BandData adapts the class template to the requirements of zero compile
// time band arguments.
*/
template<>
class BandData<>
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
   template< typename... RBAs >
   explicit inline BandData( ptrdiff_t index, RBAs... args );

   BandData( const BandData& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~BandData() = default;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   BandData& operator=( const BandData& ) = delete;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline ptrdiff_t band  () const noexcept;
   inline size_t    row   () const noexcept;
   inline size_t    column() const noexcept;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   const ptrdiff_t band_;    //!< The band index.
   const size_t    row_;     //!< The index of the row containing the first element of the band.
   const size_t    column_;  //!< The index of the column containing the first element of the band.
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief The constructor for BandData.
//
// \param index The index of the band.
// \param args The optional band arguments.
*/
template< typename... RBAs >  // Optional band arguments
inline BandData<>::BandData( ptrdiff_t index, RBAs... args )
   : band_  ( index  )                        // The band index
   , row_   ( index >= 0L ?   0UL : -index )  // The index of the row containing the first element of the band
   , column_( index >= 0L ? index :    0UL )  // The index of the column containing the first element of the band
{
   MAYBE_UNUSED( args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the index of the band of the underlying dense matrix.
//
// \return The index of the band.
*/
inline ptrdiff_t BandData<>::band() const noexcept
{
   return band_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the index of the row containing the first element of the band.
//
// \return The first row index.
*/
inline size_t BandData<>::row() const noexcept
{
   return row_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the index of the column containing the first element of the band.
//
// \return The first column index.
*/
inline size_t BandData<>::column() const noexcept
{
   return column_;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR ONE COMPILE TIME BAND INDEX
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the BandData class template for a single compile time band argument.
// \ingroup band
//
// This specialization of BandData adapts the class template to the requirements of a single
// compile time band argument.
*/
template< ptrdiff_t I >  // Compile time band index
class BandData<I>
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
   template< typename... RBAs >
   explicit inline BandData( RBAs... args );

   BandData( const BandData& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~BandData() = default;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   BandData& operator=( const BandData& ) = delete;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   static constexpr ptrdiff_t band  () noexcept;
   static constexpr size_t    row   () noexcept;
   static constexpr size_t    column() noexcept;
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief The constructor for BandData.
//
// \param args The optional band arguments.
*/
template< ptrdiff_t I >       // Compile time band index
template< typename... RBAs >  // Optional band arguments
inline BandData<I>::BandData( RBAs... args )
{
   MAYBE_UNUSED( args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the index of the band of the underlying dense matrix.
//
// \return The index of the band.
*/
template< ptrdiff_t I >  // Compile time band index
constexpr ptrdiff_t BandData<I>::band() noexcept
{
   return I;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the index of the row containing the first element of the band.
//
// \return The first row index.
*/
template< ptrdiff_t I >  // Compile time band index
constexpr size_t BandData<I>::row() noexcept
{
   return ( I >= 0L ? 0UL : -I );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the index of the column containing the first element of the band.
//
// \return The first column index.
*/
template< ptrdiff_t I >  // Compile time band index
constexpr size_t BandData<I>::column() noexcept
{
   return ( I >= 0L ? I : 0UL );
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
