///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file FeudalString.h
 * \author ghall
 * \date Oct 1, 2009
 *
 * \brief Exactly like std::string, but only 2^32 elements, uses SmallVec for containment
 */

#ifndef FEUDAL_STRING_H_
#define FEUDAL_STRING_H_

#include "Compare.h"
#include "feudal/BinaryStream.h"
#include "feudal/Mempool.h"
#include "feudal/SmallVec.h"
#include <algorithm>
#include <cstddef>
#include <iterator>
#include <limits>
#include <string>
#include <vector>

template <class charT, class Traits = std::char_traits<charT> >
class FeudalString
{
public:
    typedef Traits traits_type;
    typedef MempoolAllocator<charT> Allocator;
    typedef SmallVec<charT, Allocator> container;
    typedef charT value_type;
    typedef charT& reference;
    typedef charT const& const_reference;
    typedef charT const* const_pointer;
    typedef charT* pointer;
    typedef unsigned int size_type;
    typedef std::ptrdiff_t difference_type;
    typedef FwdIterator<charT> iterator;
    typedef FwdConstIterator<charT> const_iterator;
    typedef RevIterator<charT> reverse_iterator;
    typedef RevConstIterator<charT> const_reverse_iterator;

    static const size_type npos = -1;
    static const charT EmptyString[];

    /*
     * Constructors
     */


    /// Default constructor
    FeudalString()
    : mContainer(Allocator())
    {}

    explicit FeudalString(const Allocator& a)
    : mContainer(a)
    {}

    /// Copy constructor
    FeudalString(const FeudalString& str)
    { assign(str.begin(),str.end()); }

    /// Copy constructor: cast a std::string into a FeudalString
    template <class alloc>
    FeudalString(const std::basic_string<charT, Traits, alloc>& str,
                 Allocator const& a = Allocator())
    : mContainer(a)
    { assign(str.begin(),str.end()); }

    /// Copy constructor, copy n characters from FeudalString starting at pos.
    /// n defaults to end of string
    FeudalString(const FeudalString& str, size_type pos, size_type n = npos,
                   const Allocator& a = Allocator())
    : mContainer(a)
    { assign(str,pos,n); }

    /// Copy constructor, copy n characters from c string.
    FeudalString(const_pointer cstr, size_type n,
                    const Allocator& a=Allocator())
    : mContainer(a)
    { assign(cstr,cstr+n); }

    /// Copy constructor, copy entire null terminated c string.
    FeudalString(const_pointer cstr, const Allocator& a = Allocator())
    : mContainer(a)
    { assign(cstr,cstr+Traits::length(cstr)); }

    /// Sized constructor, create an n length string of c characters
    FeudalString(size_type n, value_type c, const Allocator& a = Allocator())
    : mContainer(a)
    { assign(n,c); }


    /// Iterator constructor, create a string from the first and last iterators
    template <class Itr>
    FeudalString(Itr first, Itr const& last, Allocator const& a = Allocator())
    : mContainer(a)
    { assign(first,last); }

    /// vector constructor, create a string from std::vector of characters
    FeudalString(const std::vector<charT>& v, Allocator const& a = Allocator())
    : mContainer(a)
    { assign(v.begin(),v.end()); }

    /// Sized constructor, create an n length string of \0 characters
    explicit FeudalString(size_type n, const Allocator& a = Allocator())
    : mContainer(a)
    { assign(n,charT()); }

    /// Sized constructor for people who are lazy about specifying constants
    /// of the correct size.
    explicit FeudalString( int n )
    { assign(n,charT()); }

    /// Single-char constructor
    explicit FeudalString( value_type c )
    { assign(1u,c); }

    /*
     * Assignment operator
     */

    /// Copy operator
    FeudalString& operator=(FeudalString const& that)
    { if ( this != &that ) assign(that.begin(),that.end()); return *this; }
    FeudalString& operator=(const_pointer cstr)
    { return assign(cstr,cstr+Traits::length(cstr)); }
    FeudalString& operator=(value_type chr)
    { return assign(1u,chr); }

    /*
     * Iterators
     */

    iterator begin() { return mContainer.begin(); }
    const_iterator begin() const { return mContainer.begin(); }
    const_iterator cbegin() const { return mContainer.cbegin(); }
    iterator begin(size_type idx) { return mContainer.begin(idx); }
    const_iterator begin(size_type idx) const { return mContainer.begin(idx); }
    const_iterator cbegin(size_type idx) const
    { return mContainer.cbegin(idx); }
    iterator end() { return mContainer.end(); }
    const_iterator end() const { return mContainer.end(); }
    const_iterator cend() const { return mContainer.cend(); }
    reverse_iterator rbegin() { return mContainer.rbegin(); }
    const_reverse_iterator rbegin() const { return mContainer.rbegin(); }
    const_reverse_iterator crbegin() const { return mContainer.crbegin(); }
    reverse_iterator rbegin(size_type idx) { return mContainer.rbegin(idx); }
    const_reverse_iterator rbegin(size_type idx) const
    { return mContainer.rbegin(idx); }
    const_reverse_iterator crbegin(size_type idx) const
    { return mContainer.crbegin(idx); }
    reverse_iterator rend() { return mContainer.rend(); }
    const_reverse_iterator rend() const { return mContainer.rend(); }
    const_reverse_iterator crend() const { return mContainer.crend(); }

    /*
     * counts and sizes
     */

    /// Return the allocator
    Allocator get_allocator() const { return mContainer.get_allocator(); }

    /// Append a character
    FeudalString& push_back(value_type val)
    { mContainer.push_back(val); return *this; }

    /// Resize the string
    FeudalString& resize(size_type sz, const_reference exemplar = charT())
    { reserve(sz); mContainer.resize(sz, exemplar); return *this; }
    /// Reserve space in the string
    FeudalString& reserve(size_type sz)
    { if ( sz ) mContainer.reserve(sz+1); return *this; }
    /// Clear the string
    FeudalString& clear() { mContainer.clear(); return *this; }
    /// Shrink the storage to fit the contents
    FeudalString& shrink_to_fit()
    { if ( mContainer.full() ) return *this;
      if ( mContainer.empty() ) mContainer.shrink_to_fit();
      else
      { mContainer.push_back(charT());
        mContainer.shrink_to_fit();
        mContainer.resize(mContainer.size()-1); }
      return *this; }

    /// Return the size of the string
    size_type size() const { return mContainer.size(); }

    size_type allocSize() const { return mContainer.size() + 1; }

    /// Return the length of the string (same as the size)
    size_type length() const { return mContainer.size(); }
    /// Return the maximum size
    size_type max_size() const { return mContainer.max_size() - 1; }
    /// Return the capacity
    size_type capacity() const { return mContainer.capacity() - 1; }
    /// Return true if the string is empty
    bool empty() const { return mContainer.empty(); }
    /// Return true if the string is at capacity
    bool full() const { return size() == capacity(); }

    /*
     * Element Access
     */

    /// Return the character at the idx position
    reference operator[](size_type idx)
    { return mContainer[idx]; }
    /// Return the character at the idx position
    const_reference operator[](size_type idx) const
    { return mContainer[idx]; }
    /// Return the character at the idx position
    reference at(size_type idx)
    { return mContainer.at(idx); }
    /// Return the character at the idx position
    const_reference at(size_type idx) const
    { return mContainer.at(idx); }

    /// Return the character at the front position
    reference front() { return mContainer.front(); }
    const_reference front() const { return mContainer.front(); }
    /// Return the character at the end position
    reference back() { return mContainer.back(); }
    const_reference back() const { return mContainer.back(); }

    /*
     * const operations
     */

    /// Return the substring of *this of length n, starting at pos
    FeudalString substr(size_type pos = 0, size_type n = npos) const
    { return FeudalString(*this,pos,n); }

    /// Return a compatible c string of *this
    const_pointer c_str() const;

    /// data
    const_pointer data() const
    { return mContainer.data(); }

    /// copy
    size_type copy(pointer dest, size_type n, size_type pos = 0) const
    { memcpy(dest,data()+pos,lenLimit(pos,n)); return n; }

    /// find
    size_type find(const_pointer cstr, size_type pos, size_type n) const;
    size_type find(value_type c, size_type pos = 0) const;
    size_type find(const FeudalString& str, size_type pos = 0) const
    { return find(str.data(), pos, str.size()); }
    size_type find(const_pointer cstr, size_type pos = 0) const
    { return find(cstr, pos, Traits::length(cstr)); }

    /// rfind
    size_type rfind(const_pointer cstr, size_type pos, size_type n) const;
    size_type rfind(value_type c, size_type pos = npos) const;
    size_type rfind(const FeudalString& str, size_type pos = npos) const
    { return rfind(str.data(), pos, str.size()); }
    size_type rfind(const_pointer cstr, size_type pos = npos) const
    { return rfind(cstr, pos, Traits::length(cstr)); }

    /// find_first_of
    size_type find_first_of(const_pointer cstr,
                                size_type pos, size_type n) const;
    size_type find_first_of(const FeudalString& str, size_type pos = 0) const
    { return find_first_of(str.data(), pos, str.size()); }
    size_type find_first_of(const_pointer cstr, size_type pos = 0) const
    { return find_first_of(cstr, pos, Traits::length(cstr)); }
    size_type find_first_of(value_type c, size_type pos = 0) const
    { return find(c, pos); }

    /// find_last_of
    size_type find_last_of(const_pointer cstr,
                                size_type pos, size_type n) const;
    size_type find_last_of(const FeudalString& str, size_type pos = npos) const
    { return find_last_of(str.data(), pos, str.size()); }
    size_type find_last_of(const_pointer cstr, size_type pos = npos) const
    { return find_last_of(cstr, pos, Traits::length(cstr)); }
    size_type find_last_of(value_type c, size_type pos = npos) const
    { return rfind(c, pos); }

    /// find_first_not_of
    size_type find_first_not_of(const_pointer cstr,
                                    size_type pos, size_type n) const;
    size_type find_first_not_of(value_type c, size_type pos = 0) const;
    size_type find_first_not_of(const FeudalString& str, size_type pos=0) const
    { return find_first_not_of(str.data(), pos, str.size()); }
    size_type find_first_not_of(const_pointer cstr, size_type pos = 0) const
    { return find_first_not_of(cstr, pos, Traits::length(cstr)); }

    /// find_last_not_of
    size_type find_last_not_of(const_pointer cstr,
                                size_type pos, size_type n) const;
    size_type find_last_not_of(value_type c, size_type pos = npos) const;
    size_type find_last_not_of(const FeudalString& str,size_type pos=npos) const
    { return find_last_not_of(str.data(), pos, str.size()); }
    size_type find_last_not_of(const_pointer cstr, size_type pos = npos) const
    { return find_last_not_of(cstr, pos, Traits::length(cstr)); }

    /// compare (modeled after the GNU STL versions)
    int compare(const FeudalString& str) const
    { using std::min;
      int result = Traits::compare(data(), str.data(), min(size(), str.size()));
      if ( !result ) result = ::compare(size(),str.size());
      return result; }
    friend int compare( FeudalString const& s1, FeudalString const& s2 )
    { return s1.compare(s2); }
    int compare(const_pointer cstr) const;
    int compare(size_type pos1, size_type n1, const FeudalString& str) const;
    int compare(size_type pos1, size_type n1, const_pointer cstr) const;
    int compare(size_type pos1, size_type n1,
                    const FeudalString& str, size_type pos2, size_type n2) const;
    int compare(size_type pos1, size_type n1,
                    const_pointer cstr, size_type n2) const;

    /*
     * modifiers
     */

    /// swap
    FeudalString& swap(FeudalString& that)
    { mContainer.swap(that.mContainer); return *this; }

    /// concatenate another string to *this
    FeudalString& operator+=(const FeudalString &str)
    { return append(str.begin(), str.end()); }
    /// concatenate another string to *this
    FeudalString& operator+=(const_pointer cstr)
    { return append(cstr, cstr + Traits::length(cstr)); }
    /// concatenate a character to *this
    FeudalString& operator+=(value_type c)
    { push_back(c); return *this; }

    /// assign
    FeudalString& assign(const FeudalString& str)
    { return assign(str.begin(),str.end()); }
    FeudalString& assign(const FeudalString& str, size_type pos, size_type n)
    { return assign(str.begin(pos),str.begin(pos+str.lenLimit(pos,n))); }
    FeudalString& assign(const_pointer cstr, size_type n)
    { return assign(cstr,cstr+n); }
    FeudalString& assign(const_pointer cstr)
    { return assign(cstr,cstr+Traits::length(cstr)); }
    FeudalString& assign(size_type n, value_type c)
    { reserve(n); mContainer.assign(n,c); return *this; }
    template <class InputIterator>
    FeudalString& assign(InputIterator first, InputIterator last)
    { using std::distance; reserve(distance(first,last));
      mContainer.assign(first,last); return *this; }

    /// append
    FeudalString& append(const FeudalString& str)
    { return append(str.begin(), str.end()); }
    FeudalString& append(const FeudalString& str, size_type pos, size_type n)
    { return append(str.begin(pos),str.begin(pos+str.lenLimit(pos,n))); }
    FeudalString& append(const_pointer cstr, size_type n)
    { return append(cstr,cstr+n); }
    FeudalString& append(const_pointer cstr)
    { return append(cstr,cstr+Traits::length(cstr)); }
    FeudalString& append(size_type sz, const value_type c)
    { reserve(size()+sz); mContainer.append(sz, c); return *this; }
    template <class Itr>
    FeudalString& append(Itr first, Itr const& last)
    { using std::distance; reserve(size()+distance(first,last));
      mContainer.append(first,last); return *this; }

    /// insert
    FeudalString& insert(size_type pos1, const FeudalString& str)
    { return insert(begin(pos1),str.begin(),str.end()); }
    FeudalString& insert(size_type pos1, const FeudalString& str,
                            size_type pos2, size_type n)
    { return insert(begin(pos1),str.begin(pos2),str.begin(pos2+str.lenLimit(pos2,n)));}
    FeudalString& insert(size_type pos1, const_pointer cstr, size_type n)
    { return insert(begin(pos1),cstr,cstr+n); }
    FeudalString& insert(size_type pos1, const_pointer cstr)
    { return insert(begin(pos1), cstr, cstr+Traits::length(cstr)); }
    FeudalString& insert(size_type pos1, size_type n, value_type c)
    { reserve(size()+n); mContainer.insert(begin(pos1),n,c); return *this; }
    iterator insert(iterator p, value_type c)
    { return mContainer.insert(p,c); }
    FeudalString& insert(iterator p, size_type n, value_type c)
    { return insert(p.pos()-data(),n,c); }
    template <class InputIterator>
    FeudalString& insert(size_type pos1,
                         InputIterator first, InputIterator last)
    { using std::distance; reserve(size()+distance(first,last));
      mContainer.insert(begin(pos1), first, last); return *this; }
    template <class InputIterator>
    FeudalString& insert(iterator p, InputIterator first, InputIterator last)
    { return insert(p.pos()-data(),first,last); }

    /// erase
    FeudalString& erase(size_type pos = 0, size_type n = npos)
    { erase(begin(pos),begin(pos+lenLimit(pos,n))); return *this; }
    iterator erase(iterator p) { return mContainer.erase(p); }
    iterator erase(iterator first, iterator last)
    { return mContainer.erase(first,last); }

    // replace
    FeudalString& replace(size_type pos1, size_type n1, const FeudalString& str)
    { return replace(begin(pos1),begin(pos1+lenLimit(pos1,n1)),
                     str.begin(),str.end()); }
    FeudalString& replace(iterator i1, iterator i2, const FeudalString& str)
    { return replace(i1,i2,str.begin(),str.end()); }
    FeudalString& replace(size_type pos1, size_type n1, const FeudalString& str,
                            size_type pos2, size_type n2)
    { return replace(begin(pos1),begin(pos1+lenLimit(pos1,n1)),
                     str.begin(pos2),str.begin(pos2+str.lenLimit(pos2,n2))); }
    FeudalString& replace(size_type pos1, size_type n1,
                                const_pointer cstr, size_type n2)
    { return replace(begin(pos1),begin(pos1+lenLimit(pos1,n1)),
                     cstr,cstr+n2); }
    FeudalString& replace(iterator i1, iterator i2,
                          const_pointer cstr, size_type n2)
    { return replace(i1,i2,cstr,cstr+n2); }
    FeudalString& replace(size_type pos1, size_type n1, const_pointer cstr)
    { return replace(begin(pos1),begin(pos1+lenLimit(pos1,n1)),
                     cstr,cstr+Traits::length(cstr)); }
    FeudalString& replace(iterator i1, iterator i2, const_pointer cstr)
    { return replace(i1,i2,cstr,cstr+Traits::length(cstr)); }
    FeudalString& replace(size_type pos1, size_type n1,
                                size_type n2, value_type c)
    { return replace(begin(pos1),begin(pos1+lenLimit(pos1,n1)),n2,c); }
    FeudalString& replace(iterator i1, iterator i2, size_type n2, value_type c)
    { for ( ; i1 != i2 && n2; ++i1, --n2 ) { *i1 = c; }
      if ( i1 != i2 ) erase(i1,i2);
      else if ( n2 ) insert(i1,n2,c);
      return *this; }
    template <class InputIterator>
    FeudalString& replace(iterator i1, iterator i2,
                          InputIterator first, InputIterator last)
    { for ( ; i1 != i2 && first != last; ++i1, ++first ) *i1 = *first;
      if ( i1 != i2 ) erase(i1,i2);
      else if ( first != last ) insert(i1,first,last);
      return *this; }


    /*
     * Broad Requirements
     */

    /// Reinitialize (aliased to clear)
    void Reinitialize() { clear().shrink_to_fit(); }
    /// Reinitialize (aliased to clear)
    void Blank() { clear().shrink_to_fit(); }

    /// Swap (aliased to swap)
    void Swap(FeudalString& that)
    { this->swap(that); }

    /// cast operator to std::string
    operator std::basic_string<charT, Traits>() const
    { return std::basic_string<charT, Traits>(c_str()); }

    /// Feudal File methods
    size_t writeBinary(BinaryWriter& writer) const
    { size_type len = size()+1;
      const_pointer ppp = empty() ? EmptyString : data();
      size_t result = writer.write(len);
      return result+writer.write(ppp,ppp+len); }

    void readBinary(BinaryReader& reader)
    { size_type len; reader.read(&len); resize(len-1);
      char buf; pointer ppp = empty() ? &buf : data();
      reader.read(ppp,ppp+len); }

    static size_t externalSizeof() { return 0; }

    size_t writeFeudal(BinaryWriter& writer, void const** ppFixed) const
    { *ppFixed = 0;
      const_pointer ppp = empty() ? EmptyString : data();
      return writer.write(ppp,ppp+size()+1); }

    /// Feudal BinaryReader constructor
    void readFeudal(BinaryReader& reader, unsigned long dataLen, void*/*fixed*/)
    { resize(dataLen-1);
      char buf; char* ppp = empty() ? &buf : data();
      reader.read(ppp,ppp+dataLen); }

    static unsigned int fixedDataLen() { return 0; }

    /// non standard STLish methods
    bool nonempty() const
    { return !mContainer.empty(); }
    void Set(const_pointer cstr, size_type i)
    { assign(cstr,i); }
    int isize() const
    { int_size_assert(); return static_cast<int>(size()); }

    /// test if a string starts with s
    bool StartsWith(const FeudalString& s) const
    { return size() >= s.size() && !Traits::compare(data(),s.data(),s.size()); }

    /// test if a string ends with s
    bool EndsWith(const FeudalString& s) const
    { return size() >= s.size() &&
        !Traits::compare(data()+size()-s.size(),s.data(),s.size()); }

    /// return true if this string can be parsed as an integer
    bool IsInt( long* = 0 ) const;
    /// return *this as an integer
    long Int() const;

    /// return true if this string can be parsed as a double
    bool IsDouble( double* = 0 ) const;
    /// return *this as a double
    double Double() const;

    /// return true if this string can be parsed as a bool
    bool IsBool( bool* = 0 ) const;
    /// return *this as a bool
    bool ToBool() const;

    /// convert *this to all upper case
    FeudalString& ToUpper();
    /// convert *this to all lower case
    FeudalString& ToLower();

    /// find a substring in *this, return -1 if not found
    int Position(const_pointer cstr) const
    { int_size_assert(); return find(cstr); }
    /// find a character in *this, return -1 if not found
    int Position(const value_type c) const
    { int_size_assert(); return find(c); }
    /// find a substring in *this, return -1 if not found
    int Position(const FeudalString& str) const
    { int_size_assert(); return find(str); }
    /// find a substring in *this up until endSearchAt, return -1 if not found
    int Position(const FeudalString& str, size_type endSearchAt) const;

    /// FirstPositionAfterRunOfAny
    size_type FirstPositionAfterRunOfAny(const_pointer cstr) const
    { return find_first_not_of(cstr); }

    /// LastPositionBeforeRunOfAny()
    int LastPositionBeforeRunOfAny(const_pointer cstr) const
    { int_size_assert(); return find_last_not_of(cstr); }

    /// PositionAfter
    int PositionAfter(const FeudalString& x, size_type startSearchAt) const
    { return find(x, startSearchAt); }

    /// return true if string contains substring at pos
    /// pos==npos means "at the end"
    bool Contains(const_pointer cstr, size_type pos) const;
    /// return true if string contains substring at pos
    /// pos==npos means "at the end"
    bool Contains(const FeudalString& str, size_type pos) const;
    /// return true if string contains substring
    bool Contains(const_pointer cstr) const
    { return find(cstr) != npos; }
    /// return true if string contains substring
    bool Contains(const FeudalString& str) const
    { return find(str) != npos; }

    /// return the part of *this before the first match for the specified string
    FeudalString Before(const_pointer cstr) const;
    /// return the part of *this before the first match for the specified string
    FeudalString Before(const FeudalString& str) const;
    /// return the part of *this before the first match for the specified string
    /// or the whole string if there is no match
    FeudalString SafeBefore(const FeudalString& str) const
    { size_type n = find(str); return (n == npos) ? *this : substr(0, n); }
    /// return the part of *this before the last match for the specified string
    /// or the whole string if there is no match
    FeudalString SafeBeforeLast(const FeudalString& str) const
    { size_type n = rfind(str); return (n == npos) ? *this : substr(0, n); }
    /// return the part of *this before the first match for the specified string
    /// or the whole string if there is no match
    FeudalString SafeBefore(const_pointer cstr) const
    { size_type n = find(cstr); return (n == npos) ? *this : substr(0, n); }

    /// return the part of *this after the first match for the specified string
    FeudalString After(const_pointer cstr) const;
    /// return the part of *this after the first match for the specified string
    FeudalString After(const FeudalString& str) const;
    /// return the part of *this after the first match for the specified string
    /// or the whole string if there is no match
    FeudalString SafeAfter(const FeudalString &str) const
    { size_type n = find(str);
      return (n == npos) ? *this : substr(n + str.size()); }
    /// return the part of *this after the last match for the specified string
    /// or the whole string if there is no match
    FeudalString SafeAfterLast(const FeudalString &str) const
    { size_type n = rfind(str);
      return (n == npos) ? *this : substr(n +str.size()); }
    /// return the part of *this after the first match for the specified string
    /// or the whole string if there is no match
    FeudalString SafeAfter(const_pointer cstr) const
    { size_type n = find(cstr);
      return (n == npos) ? *this : substr(n + Traits::length(cstr)); }

    /// return the part of *this before the final match for the specified string
    FeudalString RevBefore(const FeudalString& s) const
    { size_type pos = rfind(s); return substr(0, pos); }
    /// return the part of *this after the final match for the specified string
    FeudalString RevAfter(const FeudalString& s) const
    { size_type pos = rfind(s); return substr(pos + s.size()); }
    /// find the position of substring starting from the end of the string
    int PosRev(const FeudalString& s) const
    { int_size_assert(); return rfind(s); }

    /// return the part of *this between substring s1 and s2
    FeudalString Between(const FeudalString& s1, const FeudalString& s2) const
    { return (*this).After(s1).Before(s2); }

    /// replace the first occurrence of substring from with substring to
    FeudalString& ReplaceBy(const FeudalString& from, const FeudalString& to);

    /// replace all occurrences of substring from with substring to
    FeudalString& GlobalReplaceBy(const FeudalString& from, const FeudalString& to);

    /// return the number of occurrences of str in *this
    size_type Freq(const FeudalString& str) const;

    /// return string trimmed of all initial and final chars that occur in cstr
    FeudalString Trim(const_pointer cstr) const
    { return FeudalString(*this).TrimInPlace(cstr); }

    /// trim all initial chars that occur in cstr
    FeudalString& LTrim(const_pointer cstr)
    { erase(0,find_first_not_of(cstr)); return *this; }

    /// trim all final chars that occur in cstr
    FeudalString& RTrim(const_pointer cstr)
    { resize(find_last_not_of(cstr)+1); return *this; } // npos trick

    /// trim all initial and final chars that occur in cstr
    FeudalString& TrimInPlace(const_pointer cstr)
    { return LTrim(cstr).RTrim(cstr); }

    // ReplaceExtension()
    FeudalString ReplaceExtension(FeudalString const& ext,
                                  FeudalString const& newExt) const;

private:
    void int_size_assert() const
    { if ( size() > static_cast<size_type>(std::numeric_limits<int>::max()) )
        tooBigForInts(); }

    void tooBigForInts() const;

    size_type lenLimit(size_type pos, size_type len) const
    { AssertLe(pos,size());
      using std::min; return min(len,size()-min(size(),pos)); }

    pointer data() { return mContainer.data(); }

    // container
    container   mContainer;
};

/*
 * Operators
 */

/// operator+ (concat)
template <class charT, class Traits>
FeudalString<charT, Traits>
operator+ (const FeudalString<charT, Traits>& lhs,
           const FeudalString<charT, Traits>& rhs)
{ FeudalString<charT, Traits> result;
  result.reserve(lhs.size()+rhs.size());
  return result.append(lhs.begin(),lhs.end()).append(rhs.begin(),rhs.end()); }

template <class charT, class Traits>
FeudalString<charT, Traits>
operator+ (const charT* lhs, const FeudalString<charT, Traits>& rhs)
{ size_t len = Traits::length(lhs);
  FeudalString<charT, Traits> result;
  result.reserve(len+rhs.size());
  return result.append(lhs,lhs+len).append(rhs.begin(),rhs.end()); }

template <class charT, class Traits>
FeudalString<charT, Traits>
operator+ (const FeudalString<charT, Traits>& lhs, const charT* rhs)
{ size_t len = Traits::length(rhs);
  FeudalString<charT, Traits> result;
  result.reserve(lhs.size()+len);
  return result.append(lhs.begin(),lhs.end()).append(rhs,rhs+len); }

template <class charT, class Traits>
FeudalString<charT, Traits>
operator+ (charT lhs, const FeudalString<charT, Traits>& rhs)
{ FeudalString<charT, Traits> result;
  result.reserve(rhs.size()+1);
  return result.push_back(lhs).append(rhs.begin(),rhs.end()); }

template <class charT, class Traits>
FeudalString<charT, Traits>
operator+ (const FeudalString<charT, Traits> &lhs, charT rhs)
{ FeudalString<charT, Traits> result;
  result.reserve(lhs.size()+1);
  return result.append(lhs.begin(),lhs.end()).push_back(rhs); }

/// operator==
template <class charT, class Traits>
bool operator==(FeudalString<charT, Traits> const& v1,
                FeudalString<charT, Traits> const& v2)
{ using std::equal;
  return v1.size()==v2.size() && equal(v1.begin(), v1.end(), v2.begin()); }
template <class charT, class Traits>
bool operator==(FeudalString<charT, Traits> const& v1, charT const* cstr)
{ using std::equal;
  return v1.size()==Traits::length(cstr) && equal(v1.begin(), v1.end(), cstr); }
template <class charT, class Traits>
bool operator==(charT const* cstr, FeudalString<charT, Traits> const& v1)
{ using std::equal;
  return v1.size()==Traits::length(cstr) && equal(v1.begin(), v1.end(), cstr); }

/// operator!=
template <class charT, class Traits>
bool operator!=(FeudalString<charT, Traits> const& v1,
                FeudalString<charT, Traits> const& v2)
{ using std::equal;
  return v1.size()!=v2.size() || !equal(v1.begin(),v1.end(),v2.begin()); }
template <class charT, class Traits>
bool operator!=(FeudalString<charT, Traits> const& v1, charT const* cstr)
{ using std::equal;
  return v1.size()!=Traits::length(cstr) || !equal(v1.begin(),v1.end(),cstr); }
template <class charT, class Traits>
bool operator!=(charT const* cstr, FeudalString<charT, Traits> const& v1)
{ using std::equal;
  return v1.size()!=Traits::length(cstr) || !equal(v1.begin(),v1.end(),cstr); }

/// operator <
template <class charT, class Traits>
bool operator<(FeudalString<charT, Traits> const& v1,
               FeudalString<charT, Traits> const& v2)
{ using std::lexicographical_compare;
  return lexicographical_compare(v1.begin(),v1.end(),v2.begin(),v2.end()); }
template <class charT, class Traits>
bool operator<(FeudalString<charT, Traits> const& v1, charT const* cstr)
{ using std::lexicographical_compare;
  return lexicographical_compare(v1.begin(),v1.end(),
                                 cstr,cstr+Traits::length(cstr)); }
template <class charT, class Traits>
bool operator<(charT const* cstr, FeudalString<charT, Traits> const& v1)
{ using std::lexicographical_compare;
  return lexicographical_compare(cstr,cstr+Traits::length(cstr),
                                 v1.begin(),v1.end()); }

/// operator >
template <class charT, class Traits>
bool operator>(FeudalString<charT, Traits> const& v1,
               FeudalString<charT, Traits> const& v2)
{ return (v2 < v1); }
template <class charT, class Traits>
bool operator>(FeudalString<charT, Traits> const& v1, charT const* cstr)
{ return (cstr < v1); }
template <class charT, class Traits>
bool operator>(charT const* cstr, FeudalString<charT, Traits> const& v2)
{ return (v2 < cstr); }

/// operator <=
template <class charT, class Traits>
bool operator<=(FeudalString<charT, Traits> const& v1,
                FeudalString<charT, Traits> const& v2)
{ return !(v2 < v1); }
template <class charT, class Traits>
bool operator<=(FeudalString<charT, Traits> const& v1, charT const* cstr)
{ return !(cstr < v1); }
template <class charT, class Traits>
bool operator<=(charT const* cstr, FeudalString<charT, Traits> const& v2)
{ return !(v2 < cstr); }

/// operator >=
template <class charT, class Traits>
bool operator>=(FeudalString<charT, Traits> const& v1,
                FeudalString<charT, Traits> const& v2)
{ return !(v1 < v2); }
template <class charT, class Traits>
bool operator>=(FeudalString<charT, Traits> const& v1, charT const* cstr)
{ return !(v1 < cstr); }
template <class charT, class Traits>
bool operator>=(charT const *cstr, FeudalString<charT, Traits> const& v1)
{ return !(cstr < v1); }

template <class charT, class Traits>
void swap(FeudalString<charT, Traits>& v1, FeudalString<charT, Traits>& v2)
{ v1.swap(v2); }

template <class T, class Tr>
struct Serializability< FeudalString<T,Tr> > : public SelfSerializable {};

#endif /* FEUDAL_STRING_H_ */
