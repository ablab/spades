///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file MasterVec.h
 * \author tsharpe
 * \date Jul 28, 2009
 *
 * \brief
 */

#ifndef FEUDAL_MASTERVEC_H_
#define FEUDAL_MASTERVEC_H_


#include "feudal/IncrementalWriter.h"
#include "feudal/FeudalFileReader.h"
#include "feudal/VirtualMasterVec.h"
#include "feudal/OuterVec.h"
#include "feudal/FeudalControlBlock.h"
#include "feudal/BinaryStreamTraits.h"
#include "feudal/FeudalTools.h"
#include "math/Functions.h"
#include "random/Random.h"
#include "system/Assert.h"
#include "system/SortInPlace.h"
#include "String.h"
#include "Vec.h"
#include <cstddef>

/*
 * Additional requirements on T, the feudal type, in addition to those specified
 * by OuterVec:
 * Must have void readFeudal( BinaryReader& reader, unsigned long dataLen, void* fixed );
 * Must have size_t writeFeudal( BinaryWriter&, void const** ) const;
 * Must have static unsigned int fixedDataLen();
 */
template <class T>
class MasterVec : public OuterVec<T>
{
public:
    typedef OuterVec<T> BaseT;
    typedef typename BaseT::size_type size_type;

    MasterVec() {}

    explicit MasterVec( size_type size )
    : BaseT(size) {}

    explicit MasterVec( size_type size, T const& exemplar,
                         size_type capacity = 0 )
    : BaseT(size,exemplar,capacity) {}

    template <class Itr>
    MasterVec( Itr first, Itr const& last, size_type capacity = 0 )
    : BaseT(first,last,capacity) {}

    explicit MasterVec( String const& fileName )
    { ReadAll(fileName,true); }

    template <class S1, class A1>
    explicit MasterVec( OuterVec<T,S1,A1> const& that )
    : BaseT(that) {}

    template <class S1, class A1>
    MasterVec& operator=( OuterVec<T,S1,A1> const& that )
    { BaseT::operator=(that); return *this; }

    // compiler-supplied destructor and copying are OK

    /// This function is deprecated: Use reserve() instead.
    /// The pool size argument is ignored, anyway.
    MasterVec& Reserve( unsigned long /*raw_mem_size_ignored*/, size_type cap )
    { reserve(cap); return *this; }

    /// This function is deprecated:  Use clear().shrink_to_fit().
    MasterVec& destroy() { BaseT::clear(); BaseT::shrink_to_fit(); return *this; }

    /// This function is deprecated:  Use push_back().
    MasterVec& push_back_external( T const& val )
    { push_back(val); return *this; }

    /// This function is deprecated:  Use push_back().
    MasterVec& push_back_reserve( T const& val,
                                  size_type growthIncr = 0,
                                  float growthFact = 1.3f )
    { push_back(val,growthFact,growthIncr); return *this; }

    /// This function is deprecated:  Use append().
    MasterVec& Append( MasterVec const& that )
    { BaseT::append(that.begin(),that.end()); return *this; }

    /// This function is deprecated:  Use append().
    MasterVec& Append( MasterVec const& that, size_type from, size_type to )
    { append(that.begin(from),that.begin(to)); return *this; }

    MasterVec const& WriteAll( String const& fileName ) const
    { return WriteRange(fileName,0UL,BaseT::size()); return *this; }

    MasterVec const& WriteRange( String const& fileName,
                           size_type from, size_type to ) const
    { IncrementalWriter<T> writer(fileName.c_str(),BaseT::size());
      writer.add(BaseT::begin(from),BaseT::begin(to));
      writer.close();
      return *this; }

    MasterVec const& WriteOne( String const& fileName,
                           size_type from ) const
    { IncrementalWriter<T> writer(fileName.c_str(),BaseT::size());
      writer.add(BaseT::begin(from),BaseT::begin(from+1));
      writer.close();
      return *this; }

    /// This function replaces or appends all contents with the elements from a
    /// specified feudal file.
    MasterVec& ReadAll( String const& fileName, bool append = false )
    { if ( !append ) BaseT::clear();
      FeudalFileReader rdr(fileName.c_str());
      appendFromFeudal(rdr,0,rdr.getNElements());
      return *this; }

    /// This function appends a range of elements from a specified feudal file.
    MasterVec& ReadRange( String const& fileName,
                          size_type from, size_type to,
                          int /*extra*/ = 0 )
    { FeudalFileReader rdr(fileName.c_str());
      appendFromFeudal(rdr,from,to);
      return *this; }

    /// This function appends a specified element from a specified feudal file.
    /// It's rare that you actually need a single element put into a mastervec.
    /// Try to rework so that you use a VirtualMasterVec to directly retrieve
    /// the element you're interested in.
    MasterVec& ReadOne( String const& fileName, size_type idx )
    { FeudalFileReader rdr(fileName.c_str());
      appendFromFeudal(rdr,idx,idx+1);
      return *this; }

    /// This function appends a specified set of elements from a specified
    /// feudal file.  Note that duplicated elements are allowed and will be
    /// read as such.
    /// It's rare that you actually need a sparse, random assortment of entries.
    /// Use a VirtualMasterVec to pick off the entries you need on the fly.
    // TODO: Potentially dangerous truncation of IDs
    MasterVec& Read( String const& fileName,
                     vec<int> const& entries,
                     int /*extra*/ = 0, Bool /*pre_reserved*/ = False )
    { FeudalFileReader rdr(fileName.c_str());
      BaseT::reserve(BaseT::size()+entries.size());
      vec<int>::const_iterator end(entries.end());
      for ( vec<int>::const_iterator itr(entries.begin()); itr != end; ++itr )
      { appendFromFeudal(rdr,*itr,*itr+1); }
      return *this; }

    MasterVec& Read( String const& fileName,
                     vec<int64_t> const& entries,
                     int /*extra*/ = 0, Bool /*pre_reserved*/ = False )
    { FeudalFileReader rdr(fileName.c_str());
      BaseT::reserve(BaseT::size()+entries.size());
      typedef vec<int64_t>::const_iterator Itr;
      Itr end(entries.end());
      for ( Itr itr(entries.begin()); itr != end; ++itr )
      { appendFromFeudal(rdr,*itr,*itr+1); }
      return *this; }

    /// This function populates a specified set of elements from a specified
    /// feudal file, without disturbing any other elements, except that the
    /// MasterVec is resized to the total number of elements in the file.
    /// Now that we have good incremental reading and writing facilities you
    /// should find little need for this.
    // TODO: Potentially dangerous truncation of IDs
    MasterVec& SparseRead( String const& fileName,
                           vec<int> const& entries,
                           int /*extra*/ = 0, Bool /*pre_reserved*/ = False )
    { VirtualMasterVec<T> vVec(fileName.c_str());
      if ( BaseT::size() < vVec.size() ) resize(vVec.size());
      vec<int>::const_iterator end(entries.end());
      for ( vec<int>::const_iterator itr(entries.begin()); itr != end; ++itr )
      { BaseT::operator[](*itr) = vVec[*itr]; }
      return *this; }

    MasterVec& SparseRead( String const& fileName,
                           vec<size_t> const& entries,
                           int /*extra*/ = 0, Bool /*pre_reserved*/ = False )
    { VirtualMasterVec<T> vVec(fileName.c_str());
      if ( BaseT::size() < vVec.size() ) resize(vVec.size());
      vec<size_t>::const_iterator end(entries.end());
      for ( vec<size_t>::const_iterator itr(entries.begin()); itr != end; ++itr )
      { BaseT::operator[](*itr) = vVec[*itr]; }
      return *this; }


    MasterVec& SparseReadRange( const String& fileName,
                                size_t from, size_t to,
                                int /* extra */ = 0 )
    { VirtualMasterVec<T> vVec(fileName.c_str());
      AssertLe(to,vVec.size());
      AssertLe(from,to);
      if ( BaseT::size() < vVec.size() ) resize(vVec.size());
      for ( ; from < to; ++from ) (*this)[from] = vVec[from];
      return *this; }

    template<class BOOL_T>
    MasterVec& EraseIf( vec<BOOL_T> const& to_remove )
    { BaseT::eraseIf(to_remove.begin(),to_remove.end()); return *this; }

  
    template<class BOOL_T>
    MasterVec& EraseUnless( vec<BOOL_T> const& to_keep )
    { BaseT::eraseUnless(to_keep.begin(),to_keep.end()); return *this; }

    /// Deprecated:  Use eraseEntries.
    // TODO: Potentially dangerous truncation of IDs
    MasterVec& RemoveByIndex( vec<int> const& to_remove )
    { BaseT::eraseEntries(to_remove.begin(),to_remove.end()); return *this; }

    /// This function is deprecated.  Just use swap.
    void SwapElements( size_type iii, size_type jjj )
    { BaseT::operator[](iii).swap(BaseT::operator[](jjj)); }

    /// It's rare you'd need this.  We don't need to pre-allocate raw
    /// data pools anymore.
    static size_type MastervecFileRawCount( const String& fileName )
    { FeudalControlBlock fcb(fileName.c_str());
      return fcb.getVarDataLen()/sizeof(typename T::value_type); }

    // ports from FeudalTemplate
    template<class V>
    void ElementSizes(V& sizelist) const
    { typedef typename BaseT::const_iterator Itr;
      sizelist.clear();
      sizelist.reserve(BaseT::size());
      for ( Itr itr(BaseT::begin()), end(BaseT::end()); itr != end; ++itr )
        sizelist.push_back(itr->size()); }

    size_t SizeSum() const
    { size_t result = 0;
      typedef typename BaseT::const_iterator Itr;
      for ( Itr itr(BaseT::begin()), end(BaseT::end()); itr != end; ++itr )
          result += itr->size();
      return result; }

    // Return the maximum size of an entry.

    size_t MaxSize( ) const
    { typedef typename BaseT::const_iterator Itr;
      using std::max;
      typename T::size_type result = 0;
      for ( Itr itr(BaseT::begin()), end(BaseT::end()); itr != end; ++itr )
          result =  max(result,itr->size());
      return result; }

    // given a list v of positive integers, suppose that each entry n is
    // replaced by n copies of itself.  Then the median of the resulting
    // enlarged list is the "N50" of the original list.
    int N50() const
    { vec<typename T::size_type> eltSizes;
      ElementSizes(eltSizes);
      ::Sort(eltSizes);
      return ::N50(eltSizes); }

    // Sort this vector in place.
    void Sort() { sortInPlace(this->begin(),this->end()); }

    // Sort this vector in place.  perm is permuted to the same sort order.
    void SortSync( std::vector<size_type> & perm )
    { if ( BaseT::size() > 1 )
      { QuickSort(0, this->size() - 1, perm); InsertionSort(perm); } }

    // Cull all adjacent duplicates.  Done on a sorted array, this
    // will result in only unique elements in the entire set.
    void Unique()
    { if (this->empty()) return; // EARLY RETURN!
      typedef typename BaseT::iterator Itr;
      using std::iter_swap;
      Itr dest(BaseT::begin());
      for ( Itr itr(dest+1), end(BaseT::end()); itr != end; ++itr )
          if ( !(*itr == *dest) ) iter_swap(itr,++dest);
      BaseT::resize(dest-BaseT::begin()+1); }

    // Sort this list, and remove all duplicate items.
    void UniqueSort() { Sort(); Unique(); }

    /// Returns a very large number -- we don't really have a fixed raw
    /// capacity any more.  This can go away when mastervec goes away, and
    /// FeudalDataManagerTemplate.h can be adjusted accordingly.  It's the
    /// only thing that refers to this useless method.
    size_t rawcapacity() const { return ~0UL >> 1; }

private:
    void preAlloc( FeudalFileReader const& rdr, size_type start, size_type end )
    { typename T::value_type* ppp = static_cast<typename T::value_type*>(0);
      size_t extSize = BinaryReader::externalSizeof(ppp);
      if ( extSize )
      { size_t maxSize = BaseT::getSubAllocator().getMaxEnchunkableSize();
        size_t preAllocSize = rdr.getPreallocation(maxSize,extSize,start,end);
        BaseT::getSubAllocator().preAllocate(1UL,preAllocSize); } }

    void appendFromFeudal( FeudalFileReader& rdr, size_t start, size_t stop )
    {
        AssertLe(start,stop);
        AssertLe(stop,rdr.getNElements());
        size_t dist = stop - start;
        BaseT::resize(BaseT::size()+dist);
        preAlloc(rdr,start,stop);
        typedef typename BaseT::iterator Itr;
        Itr end(BaseT::end());
        for ( Itr itr(end-dist); itr != end; ++itr, ++start )
            itr->readFeudal(rdr.getData(start),rdr.getDataLen(start),
                            rdr.getFixedData(start,T::fixedDataLen()));
    }

    // quicksort, will only sort chunks larger than qsort_cutoff.
    // InsertionSort takes care of the rest.  perm is sorted with
    // the same sort order as this.
    void QuickSort(size_type lower, const size_type upper,
            std::vector<size_type>& perm)
    {   using std::swap;
        if ( upper < lower + QSORT_CUTOFF)
            return;
        do
        {
            size_type i = lower;
            size_type j = upper + 1;
            size_type k = lower + randomx()%(j-i);
            (*this)[i].swap((*this)[k]);
            swap(perm[i], perm[k]);
            while ( 1 )
            {   do ++i; while ( i <= upper && (*this)[i] < (*this)[lower] );
                do --j; while ( (*this)[lower] < (*this)[j] );
                if (i > j) break;
                (*this)[i].swap((*this)[j]);
                swap(perm[i], perm[j]); }
            (*this)[lower].swap((*this)[j]);
            swap(perm[lower], perm[j]);
            if ( j > lower + QSORT_CUTOFF ) QuickSort(lower, j-1, perm);
            lower = j + 1;
        }
        while ( upper >= lower + QSORT_CUTOFF ); }

    // InsertionSort, more efficient at sorting smaller chunks in a
    // nearly-sorted array. Perm is sorted with the same sort order as this.
    void InsertionSort(std::vector<size_type> & perm)
    {   using std::swap;
        for (size_type i = 1; i < this->size(); ++i)
            for (size_type j = i; j > 0 && (*this)[j] < (*this)[j-1]; --j)
            {   (*this)[j-1].swap((*this)[j]);
                swap(perm[j-1], perm[j]); } }

    static const size_type QSORT_CUTOFF = 20;
};

template <class T>
struct Serializability< MasterVec<T> > : public SelfSerializable {};

/// Deprecated:  Use v.clear().shrink_to_fit()
template <class T>
inline void Destroy( MasterVec<T>& v )
{ v.destroy(); }

// Return the position of an element in a sorted mastervec, else -1.  If the
// element appears more than once, the position of one of its instances is
// returned.
template<class T>
int BinPosition( const MasterVec<T>& mvec, const T& item)
{
    if (mvec.size() == 0)
        return -1;

    size_type first = 0;
    size_type last = mvec.size() - 1;
    size_type next;
    while (1)
    {
        if (first == last)
            return (!(item < mvec[last]) && !(mvec[last] < item)) ? last : -1;
        next = first + (last - first) / 2;
        if (item < mvec[next])
            last = next;
        else if (mvec[next] < item)
            first = next + 1;
        else
            return next;
    }
}

#endif /* FEUDAL_MASTERVEC_H_ */
