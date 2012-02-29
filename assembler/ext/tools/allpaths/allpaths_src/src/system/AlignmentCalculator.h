///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file AlignmentCalculator.h
 * \author tsharpe
 * \date Jun 21, 2010
 *
 * \brief Calculator for alignment requirements and padding for any class.
 */
#ifndef ALIGNMENTCALCULATOR_H_
#define ALIGNMENTCALCULATOR_H_

#include <cstddef>

template <class T>
class AlignmentCalculator
{
    class T0
    { char mC; T mT; };
    class T1 : T
    {
      public:
        size_t offset() const
        { return &mC - reinterpret_cast<char const*>(this); }
      private:
        char mC;
    };

public:
    static size_t getAlignment()
    { return sizeof(T0) - sizeof(T); }
    static size_t getTailPadding()
    { return sizeof(T)-reinterpret_cast<T1*>(0)->offset(); }
};

template <> inline size_t AlignmentCalculator<size_t>::getTailPadding()
{ return 0; }

#endif /* ALIGNMENTCALCULATOR_H_ */
