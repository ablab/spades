///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "feudal/BaseVec.h"
#include "graph/Digraph.h"
#include "graph/DigraphTemplate.h"
#include "dna/Bases.h"
#include "system/Assert.h"
#include "system/ErrNo.h"
#include "system/Exit.h"
#include "system/file/FileReader.h"
#include <algorithm>
#include <cstring>
#include <iostream>
#include <unistd.h>

using std::istream;
using std::ostream;
using std::cout;
using std::endl;

// Cap: in a given basevector, replace any sequence of N > n identical
// bases by n of the same base.
void BaseVec::Cap( unsigned nnn )
{
    if ( nnn > size() )
        return;

    iterator dest(begin());
    iterator src(begin());
    iterator stop(end());
    value_type lastBase = 0;
    unsigned count = 0;
    while ( src != stop )
    {
        value_type base = *src;
        if ( base != lastBase )
        {
            count = 0;
            lastBase = base;
        }
        if ( ++count <= nnn )
        {
            if ( src != dest )
                dest.set( base );
            ++dest;
        }
        ++src;
    }
    resize(dest.pos());
}

// Return the first position of other inside ourselves, or size() if
// not found.  Start looking at start, and stop looking when we reach
// the minimum of end or size().
BaseVec::size_type BaseVec::Find(const BaseVec& other,
                                    BaseVec::size_type start,
                                    BaseVec::size_type end) const
{
    AssertLe(end,size());
    AssertLe(start,end);

    if (other.size() > end-start) { return size(); }

    const_iterator my_end(begin(end-other.size()+1));
    const_iterator my_iter = begin(start);
    const_iterator other_end(other.end());
    for ( ; my_iter != my_end; ++my_iter )
    {
        bool found = true;
        const_iterator other_iter(other.begin());
        const_iterator my_iter2(my_iter);
        for ( ; found && other_iter != other_end; ++other_iter, ++my_iter2 )
            found = (*my_iter2 == *other_iter);

        if ( found )
            return my_iter.pos();
    }

    return size();
}

// FindAll: Find all start positions of "other" inside "this".
vec<BaseVec::size_type> BaseVec::FindAll(const BaseVec& other) const
{
    vec<size_type> places;
    size_type nnn = size();
    size_type idx = 0;
    while ( (idx = Find(other,idx,nnn)) != nnn )
    {
        places.push_back(idx++);
    }
    return places;
}

// Return true if two basevectors overlap exactly by r bases.
// i.e., the last r bases of this are equal to the first r bases of that.
bool BaseVec::Overlap(const BaseVec& that, BaseVec::size_type r) const
{
    if ((that.size() < r) || (size() < r)) { return false; }
    const_iterator stop(end());
    const_iterator itr(begin(size()-r));
    const_iterator itr2(that.begin());
    bool result = true;
    for ( ;result && itr != stop; ++itr, ++itr2 )
        result = (*itr == *itr2);
    return result;
}

// SetToSubOf(orig_bv, start_pos, len):  Set this to the length len
// sub-basevector of orig_bv, starting at position start_pos.  The
// case where this == &orig_bv is allowed.  If len is -1, it's
// adjusted to mean "the end" of the bvec being copied.
BaseVec& BaseVec::SetToSubOf(const BaseVec& source,
                BaseVec::size_type start, BaseVec::size_type len,
                BaseVec::size_type extra)
{
    AssertLe( start, source.size() );
    if ( len == ~0U )
        len = source.size() - start;
    AssertLe( len, source.size()-start );

    reserve(len+extra);
    setSize(len);

    if ( (len = (len + 3) / 4) )
    {
        value_type const* src = source.data() + start/4;
        value_type* dst = data();
        switch ( start & 3 )
        {
        case 0:
            if ( dst != src )
                memcpy(dst,src,len);
            break;
        case 1:
            while ( len-- )
            {
                *dst++ = (src[0] >> 2) | (src[1] << 6);
                src += 1;
            }
            break;
        case 2:
            while ( len-- )
            {
                *dst++ = (src[0] >> 4) | (src[1] << 4);
                src += 1;
            }
            break;
        case 3:
            while ( len-- )
            {
                *dst++ = (src[0] >> 6) | (src[1] << 2);
                src += 1;
            }
            break;
        }
    }

    return *this;
}

//
// Genomic Classifications
//

// Return true if all bases are equal or if empty.
bool BaseVec::IsHomopolymer() const
{
    bool result = true;
    if ( !empty() )
    {
        const_iterator stop(end());
        const_iterator itr(begin());
        value_type firstBase = *itr;
        for ( ++itr; result && itr != stop; ++itr )
            result = (*itr == firstBase);
    }
    return result;
}

// Return the % of highest base and which base it is.
pair<float, unsigned char> BaseVec::HomopolPercent() const
{
    pair<float, unsigned char> result(-1, 255);
    if ( !empty() )
    {
        unsigned int count[4] = {0,0,0,0};
        const_iterator stop(end());
        for ( const_iterator itr(begin()); itr != stop; ++itr )
            ++count[*itr];

        unsigned int* mx = max_element(count, count+4);
        result.first = *mx * 100.0f / size();
        result.second = mx - count;
    }
    return result;
}

// Return homopolymer count at this base (extends and counts both
// left and right from the specified base).
int BaseVec::Homopol(BaseVec::size_type idx) const
{
    AssertLe(idx,size());
    const_iterator pos(begin(idx));
    const value_type firstBase = *pos;
    int ret = 1;

    const_iterator stop(end());
    const_iterator itr(pos);
    while( ++itr != stop && *itr == firstBase )
        ret += 1;

    const_iterator start(begin());
    while ( itr != start && *--pos == firstBase )
        ret += 1;

    return ret;
}

// Replace *this by its reverse complement.
BaseVec& BaseVec::ReverseComplement()
{
    if ( !(size() & 3) ) // if size is evenly divisible by 4
    {
        // run algorithm byte-wise
        value_type* head = data();
        value_type* tail = dataEnd();
        while ( head != tail )
        {
            value_type tmp = Base::rcByte(*--tail);
            if ( head == tail )
            {
                *head = tmp;
                break;
            }
            *tail = Base::rcByte(*head);
            *head++ = tmp;
        }
    }
    else
    {
        // run algorithm base-wise
        iterator head = begin();
        iterator tail = end();
        while (head != tail)
        {
            value_type tmp = Complement(*--tail);
            if (head == tail)
            {
                head.set(tmp);
                break;
            }
            tail.set(Complement(*head));
            head.set(tmp);
            ++head;
        }
    }
    return *this;
}

// Return the CanonicalForm (without altering the bvec).
//  - if this is less than its RC, return fw_form
//  - if RC is less than this, return rc_form
//  - if this is equal to its RC, return palindrome
BaseVec::CanonicalForm BaseVec::getCanonicalForm() const
{
    const_iterator fwd(begin());
    const_iterator rev(end());
    while ( fwd <= --rev )
    {
        unsigned char f = *fwd;
        unsigned char r = Complement(*rev);
        if ( f < r )
            return fw_form;
        if ( r < f )
            return rc_form;
        ++fwd;
    }
    return palindrome;
}

BaseVec::size_type BaseVec::GcBases( BaseVec::size_type start,
            BaseVec::size_type end ) const
{
    AssertLe(start,size());
    AssertLe(end,size());

    size_type gc = 0;
    const_iterator tail = begin(end);
    for ( const_iterator head(begin(start)); head != tail; ++head )
        if( IsGC(*head) ) gc += 1;

    return gc;
}

unsigned int BaseVec::extractKmer( unsigned int offset, unsigned int k ) const
{
    AssertLe(k,16u);
    AssertLe(k,size());
    AssertLe(offset,size()-k);

    unsigned int result = 0;
    offset += k; // point to the end
    while ( k && (offset&3) ) // while we're not on a byte boundary
    {
        result = (result << 2) | (*this)[--offset];
        k -= 1;
    }
    if ( k >= 4 ) // if there's a whole byte to be gotten, load bytes
    {
        unsigned char const* pData = data() + offset/4;
        offset -= 4 * (k/4);
        do
        {
            result = (result << 8) | *--pData;
        }
        while ( (k-= 4) >= 4 );
    }
    while ( k-- )
        result = (result << 2) | (*this)[--offset];

    return result;
}

unsigned int BaseVec::hash( unsigned int byteOffset, unsigned int nBytes ) const
{
    AssertLe(nBytes,physicalSize(size()));
    AssertLe(byteOffset,physicalSize(size())-nBytes);

    // This is a 32-bit FNV-1a hash.  Google it.
    unsigned int hash = 2166136261;
    const unsigned int prime = 16777619;

    unsigned char const* itr = data() + byteOffset;
    unsigned char const* end = itr + nBytes;
    while ( itr != end )
    {
      hash ^= *itr++;
      hash *= prime;
    }
    return hash;
}

//
// Output
//

// Translate basevector to String of base letters (ACGT)
String BaseVec::ToString() const
{
    String s; s.reserve(size());
    std::transform(begin(),end(),std::back_inserter(s),BaseToCharMapper());
   return s;
}

//
// Global functions related to BaseVec
//

// Algorithm: Step through basevector s.  At each location, look for an
// overlap - i.e., a perfect match with the beginning of t, starting at that
// location in s, and continuing to the end of s.  If we find an overlap,
// it must be the largest overlap (because we've been searching from the
// beginning of s.)
unsigned int LargestOverlap(const BaseVec& s, const BaseVec& t,
                            unsigned int r_max, unsigned int r_min)
{
    AssertLe(r_max,s.size());
    AssertLe(r_min,r_max);

    if ( r_max > t.size() ) r_max = t.size();
    for ( ; r_max >= r_min; --r_max )
        if ( s.Overlap(t,r_max) )
            return r_max;
    return 0U;
}

// Copy count bases from src at src_start to target at target_start.
// if rc_form is true, the src will be RC'd before it is copied.
void CopyBases(const BaseVec& src, BaseVec::size_type src_start,
               BaseVec& target, BaseVec::size_type target_start,
               BaseVec::size_type count, Bool rc_from)
{
    AssertLe( src_start, src.size( ) );
    AssertLe( count, src.size()-src_start );
    AssertLe( target_start, target.size() );
    AssertLe( count, target.size()-target_start );

    // TODO: Byte level optimization
    BaseVec::iterator dst(target.begin(target_start));
    if ( rc_from )
    {
        BaseVec::const_rc_iterator end(src.rcbegin(src_start+count));
        BaseVec::const_rc_iterator itr(src.rcbegin(src_start));
        for ( ; itr != end; ++itr, ++dst )
            dst.set(*itr);
    }
    else
    {
        BaseVec::const_iterator end(src.begin(src_start+count));
        BaseVec::const_iterator itr(src.begin(src_start));
        for ( ; itr != end; ++itr, ++dst )
            dst.set(*itr);
    }
}

void BinaryWrite(int fd, const BaseVec& b)
{
    BaseVec::size_type len = b.size();
    if ( write(fd,&len,sizeof(len)) != sizeof(len) )
    {
        ErrNo err;
        cout << "BinaryWrite of Basevector length failed" << err << endl;
        CRD::exit(1);
    }
    BinaryWriteContent(fd,b);
}

void BinaryWriteContent( int fd, BaseVec const& b )
{
    BaseVec::size_type physLen = b.physicalSize(b.size());
    if ( write(fd,b.data(),physLen) != physLen )
    {
        ErrNo err;
        cout << "BinaryWrite of Basevector data failed. " << err << endl;
        CRD::exit(1);
    }
}
void BinaryRead(int fd, BaseVec& b)
{
    BaseVec::size_type len;
    if ( read(fd,&len,sizeof(len)) != sizeof(len) )
    {
        ErrNo err;
        cout << "BinaryRead of Basevector length failed. " << err << endl;
        CRD::exit(1);
    }
    b.resize(len);
    BinaryReadContent(fd,b);
}

void BinaryReadContent( int fd, BaseVec& b )
{
    BaseVec::size_type physLen = b.physicalSize(b.size());
    if ( read(fd,b.data(),physLen) != physLen )
    {
        ErrNo err;
        cout << "BinaryRead of Basevector data failed. " << err << endl;
        CRD::exit(1);
    }
}

bool BaseVec::IsGoodFeudalFile(const String& filename, bool verbose) {
    unsigned long fileSize;
    FeudalControlBlock fcb(filename.c_str(),false,&fileSize);
    bool result = true;
    if ( !fcb.isValid(filename.c_str(),fileSize,verbose) )
        result = false;
    else if ( fcb.getSizeofFixed() &&
              fcb.getSizeofFixed() != BaseVec::fixedDataLen() )
    {
        result = false;
        if ( verbose )
            cout << "Feudal file for basevector has a fixed-data length of " <<
                    fcb.getSizeofFixed() << " but we expect " <<
                    BaseVec::fixedDataLen() << "." << endl;
    }
    else if ( fcb.getSizeofX() &&
              fcb.getSizeofX()!=12 &&
              fcb.getSizeofX()!=sizeof(BaseVec) )
    {
        result = false;
        if ( verbose )
            cout << "Feudal file for basevector has a vector size of " <<
                    fcb.getSizeofX() << " but we expect 24." << endl;
    }
    else if ( fcb.getSizeofA() &&
              fcb.getSizeofA()!=4 &&
              fcb.getSizeofA()!=sizeof(BaseVec::value_type) )
    {
        result = false;
        if ( verbose )
            cout << "Feudal file for basevector has a value size of " <<
                    fcb.getSizeofA() << " but we expect " <<
                    sizeof(BaseVec::value_type) << "." << endl;
    }
    return result;
}

std::vector<bvec::size_type> BaseVec::getSizes( String const& fastbFilename )
{
    if ( !IsGoodFeudalFile(fastbFilename,true) )
        FatalErr("Feudal file " << fastbFilename << " looks fishy.");

    std::vector<bvec::size_type> results;
    FeudalFileReader ffr(fastbFilename.c_str());
    size_t nnn = ffr.getNElements();
    if ( nnn )
    {
        results.resize(nnn);
        void* fixedTable = ffr.getFixedData(0,sizeof(bvec::size_type));
        memcpy(&results.front(),fixedTable,nnn*sizeof(bvec::size_type));
    }
    return results;
}

#include "feudal/FieldVecDefs.h"
template class FieldVec< 2, MempoolAllocator<unsigned char> >;
template void digraphE<BaseVec>::AddEdge(int, int, BaseVec const&);
template void digraphE<BaseVec>::DeleteEdgeFrom(int, int);
template void digraphE<BaseVec>::AddVertices(int);
template void digraphE<BaseVec>::DeleteEdges(vec<int> const&);
template int digraphE<BaseVec>::EdgeObjectIndexByIndexTo(int, int) const;
template const BaseVec& digraphE<BaseVec>::EdgeObject(int) const;
template void digraphE<BaseVec>::ComponentsE(vec<vec<int> >&) const;
