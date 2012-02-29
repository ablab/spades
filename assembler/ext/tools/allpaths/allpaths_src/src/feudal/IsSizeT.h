///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file IsSizeT.h
 * \author tsharpe
 * \date Oct 19, 2009
 *
 * \brief
 */
#ifndef ISSIZET_H_
#define ISSIZET_H_

struct YesSizeT {};
struct NoSizeT {};

template <class T>
struct IsSizeT : public NoSizeT {};

template<>
struct IsSizeT<int> : public YesSizeT {};

template<>
struct IsSizeT<unsigned int> : public YesSizeT {};

template<>
struct IsSizeT<long> : public YesSizeT {};

template<>
struct IsSizeT<unsigned long> : public YesSizeT {};

#endif /* ISSIZET_H_ */
