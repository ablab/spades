///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file Charvector.h
 * \author jbutler
 * \date June 25, 2004
 *
 * \brief Feudal vectors of chars, uchars, and Bools.
 */
#ifndef CHARVECTOR_H_
#define CHARVECTOR_H_

#include "feudal/MasterVec.h"
#include "feudal/SerfVec.h"

typedef SerfVec<char> charvector;
typedef MasterVec< charvector > veccharvector;

typedef SerfVec<unsigned char> ucharvector;
typedef MasterVec< ucharvector > vecucharvector;

typedef ucharvector Boolvector;
typedef vecucharvector vecBoolvector;

// Remove newlines from in and save the result in out.
void StripNewlines( const charvector &in, charvector &out );

#endif /* CHARVECTOR_H_ */
