///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file IncrementalWriter.h
 * \author tsharpe
 * \date Aug 18, 2009
 *
 * \brief
 */

#ifndef FEUDAL_INCREMENTALWRITER_H_
#define FEUDAL_INCREMENTALWRITER_H_

#include "feudal/FeudalFileWriter.h"
#include <cstddef>
#include <iterator>

/// A class for writing feudal files for some particular feudal type.
/// The template parameter T must be a class that has a value_type, the
/// method size_t T::writeFeudal( BinaryWriter&, void const** ) const, and the
/// static method T::fixedDataLen().
template <class T>
class IncrementalWriter
{
public:
    IncrementalWriter( char const* filename,
                            unsigned long estimatedSize = 1000000 )
    : mWriter(filename,sizeof(T),sizeof(typename T::value_type),
              T::fixedDataLen(),estimatedSize)
    {}

    /// Construct from something with a c_str() member (like string or String)
    template <class C>
    explicit IncrementalWriter( C const& filename,
                                    unsigned long estimatedSize = 1000000,
                                    char const*(C::*)() const=&C::c_str )
    : mWriter(filename.c_str(),sizeof(T),sizeof(typename T::value_type),
              T::fixedDataLen(),estimatedSize)
    {}

    // compiler-supplied destructor is OK

    void add( T const& val )
    { size_t buf;
      void const* pBuf = &buf;
      size_t len = val.writeFeudal(mWriter.getWriter(),&pBuf);
      mWriter.addElement(len,pBuf); }

    template <class Itr>
    void add( Itr begin, Itr end )
    { while ( begin != end ) { add(*begin); ++begin; } }

    void checkPoint() { mWriter.checkPoint(); }

    unsigned long getNElements() const { return mWriter.getNElements(); }

    void close() { mWriter.close(); }

    class iterator : public std::iterator<std::output_iterator_tag,T>
    {
    public:
        iterator() : mpWriter(0) {}
        iterator( IncrementalWriter& writer ) : mpWriter(&writer) {}

        // compiler-supplied copying and destructor are OK

        iterator& operator=( T const& val )
        { mpWriter->add(val); return *this; }

        iterator& operator*() { return *this; }
        iterator& operator++() { return *this; }
        iterator& operator++(int) { return *this; }

        friend bool operator==( iterator const& itr1, iterator const& itr2 )
        { return itr1.mpWriter == itr2.mpWriter; }
        friend bool operator!=( iterator const& itr1, iterator const& itr2 )
        { return !(itr1 == itr2); }

    private:
        IncrementalWriter* mpWriter;
    };

    iterator begin() { return iterator(*this); }
    iterator end() { return iterator(); }

private:
    IncrementalWriter( IncrementalWriter const& ); // unimplemented -- no copying
    IncrementalWriter& operator=( IncrementalWriter const& ); // unimplemented -- no copying

    FeudalFileWriter mWriter;
};

#endif /* FEUDAL_INCREMENTALWRITER_H_ */
