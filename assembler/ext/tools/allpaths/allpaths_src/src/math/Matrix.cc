///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "system/Types.h"
#include "math/Matrix.h"
#include "math/MatrixTemplate.h"

MATRIX_DEF(float)
MATRIX_DEF(double)
MATRIX_DEF(Bool)

template matrix<char>::matrix( int r, int c );
template matrix<char>::matrix( int r, int c, const char& t );

template< > void matrix<char>::Print( ostream& out ) const
{    for ( int i = 0; i < Nrows( ); i++ )
     {    for ( int j = 0; j < Ncols( ); j++ )
               out << (*this)(i, j);
          out << "\n";    }    }
