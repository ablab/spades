////////////////////////////////////////////////////////////////////////////
//                  SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//      This software and its documentation are copyright (2009) by the   //
//  Broad Institute.  All rights are reserved.  This software is supplied //
//  without any warranty or guaranteed support whatsoever. The Broad      //
//  Institute is not responsible for its use, misuse, or functionality.   //
////////////////////////////////////////////////////////////////////////////
/*
 * \file QualNibbleVec.cc
 * \author tsharpe
 * \date Oct 21, 2009
 *
 * \brief
 */
#include "feudal/QualNibbleVec.h"
#include "feudal/FieldVecDefs.h"
#include "feudal/OuterVecDefs.h"

template class FieldVec< 4, MempoolAllocator<unsigned char> >;
template class OuterVec<QualNibbleVec>;

QualNibbleVec& QualNibbleVec::squash( unsigned radius )
{
    size_t nnn = size();
    value_type* vals = new value_type[nnn];
    value_type* end = vals + nnn;

    value_type* ppp = vals;
    for ( size_t idx = 0; idx < nnn; ++idx )
        *ppp++ = operator[](idx);

    using std::max;
    using std::min;
    using std::min_element;
    ppp = vals;
    for ( size_t idx = 0; idx < nnn; ++idx, ++ppp )
        set(idx,*min_element(max(vals,ppp-radius),min(ppp+radius+1,end)));

    delete [] vals;
    return *this;
}

void WriteAll(QualNibbleVecVec const& quals, String const& fn)
{
  IncrementalWriter<QualVec> writer(fn.c_str());
  for (size_t i = 0; i < quals.size(); ++i)
    writer.add(quals[i].GetQualvector());
  writer.close();
}

void LoadQualNibbleVec(const String & fn,
                       QualNibbleVecVec * quals)
{
  typedef VirtualMasterVec<QualVec> VVQV;
  typedef VVQV::const_iterator Itr;

  VVQV vvqv(fn.c_str());
  quals->clear().reserve(vvqv.size());
  for (Itr itr(vvqv.begin()), end(vvqv.end()); itr != end; ++itr )
    quals->push_back(QualNibbleVec(*itr));
}
