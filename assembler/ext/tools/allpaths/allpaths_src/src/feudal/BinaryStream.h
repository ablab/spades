///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file BinaryStream.h
 * \author tsharpe
 * \date Aug 12, 2009
 *
 * \brief Facility for serializing data structures to a file.
 * This is an incredibly naive version that doesn't handle any of the niceties
 * such as object cycles, version compatibility, etc.
 */
#ifndef FEUDAL_BINARYSTREAM_H_
#define FEUDAL_BINARYSTREAM_H_

#include "feudal/BinaryStreamTraits.h"
#include "system/Assert.h"
#include "system/file/FileReader.h"
#include "system/file/FileWriter.h"
#include <cstddef>
#include <cstring>
#include <iterator>
#include <string>
#include <utility>
#include <vector>

class MagicToken
{
public:
    MagicToken()
    { memcpy(mToken,"BINWRITE",sizeof(mToken)); }

    // compiler-supplied destructor and copying are OK

    bool isValid()
    { return !memcmp(mToken,"BINWRITE",sizeof(mToken)); }

private:
    char mToken[8];
};
TRIVIALLY_SERIALIZABLE(MagicToken);

/// Writer of binary streams.
class BinaryWriter
{
public:
    BinaryWriter( char const* filename, bool writeHeader=true )
    : mFW(filename), mpBuf(mBuf)
    { if ( writeHeader ) write(MagicToken()); }

    /// Construct from something with a c_str() member (like string or String)
    template <class C>
    explicit BinaryWriter( C const& filename, bool writeHeader=true,
                                char const*(C::*)() const=&C::c_str )
    : mFW(filename.c_str()), mpBuf(mBuf)
    { if ( writeHeader ) write(MagicToken()); }

    BinaryWriter( int fd, char const* pseudoFilename )
    : mFW(fd,pseudoFilename), mpBuf(mBuf)
    {}

    ~BinaryWriter()
    { close(); }

    std::string const& getFilename() const { return mFW.getFilename(); }

    /// Write something.
    template <class T>
    size_t write( T const& val )
    { return writeVal(val,Serializability<T>()); }

    /// Write an array of somethings.
    template <class T>
    size_t write( T const* begin, T const* end )
    { return writeArray(begin,end,Serializability<T>()); }

    /// Write a sequence from an iterator pair.
    template <class Itr>
    size_t writeItr( Itr begin, Itr const& end )
    { size_t len = 0;
      while ( begin != end ) { len += write(*begin); ++begin; }
      return len; }

    void close()
    { flush(); mFW.close(); }

    bool isOpen() const { return mFW.isOpen(); }

    void seek( size_t pos )
    { flush(); mFW.seek(pos); }

    void flush()
    { if ( mpBuf != mBuf ) { write(mBuf,mpBuf-mBuf); mpBuf = mBuf; } }

    /// Write a single object to a file in a single, simple function call.
    template <class T>
    static void writeFile( char const* filename, T const& obj )
    { BinaryWriter writer(filename); writer.write(obj); writer.close(); }

    /// Write a single object to a file in a single, simple function call.
    /// The file is described by something with a c_str() member.
    template <class C, class T>
    static void writeFile( C const& filename, T const& obj,
                                char const*(C::*)() const=&C::c_str )
    { writeFile(filename.c_str(),obj); }

protected:

    template <class SIZ_T>
    void patchSize( SIZ_T val )
    { flush(); seek(sizeof(MagicToken)); writeBuf(&val,sizeof(val)); }

private:
    BinaryWriter( BinaryWriter const& ); // unimplemented -- no copying
    BinaryWriter& operator=( BinaryWriter const& ); // unimplemented -- no copying

    template <class T>
    size_t writeVal( T const& val, TriviallySerializable )
    { writeBuf(&val,sizeof(T)); return sizeof(T); }

    template <class T>
    size_t writeVal( T const& val, SelfSerializable )
    { return val.writeBinary(*this); }

    template <class T>
    size_t writeVal( T const& val, ExternallySerializable )
    { return writeBinary(*this,val); }

    template <class T>
    size_t writeArray( T const* begin, T const* end, TriviallySerializable )
    { size_t len = sizeof(T)*(end-begin); writeBuf(begin,len); return len; }

    template <class T>
    size_t writeArray( T const* begin, T const* end, SelfSerializable )
    { size_t len = 0;
      while ( begin != end ) len += write(*begin++);
      return len; }

    template <class T>
    size_t writeArray( T const* begin, T const* end, ExternallySerializable )
    { size_t len = 0;
      while ( begin != end ) len += write(*begin++);
      return len; }

    void writeBuf( void const* data, size_t len )
    { size_t remain = mBuf+sizeof(mBuf)-mpBuf;
      if ( remain > len )
      { memcpy(mpBuf,data,len); mpBuf += len; }
      else
      { memcpy(mpBuf,data,remain);
        write(mBuf,sizeof(mBuf));
        mpBuf = mBuf;
        if ( (len -= remain) )
        { data = static_cast<char const*>(data)+remain;
          if ( len >= sizeof(mBuf) ) write(data,len);
          else { memcpy(mBuf,data,len); mpBuf += len; } } } }

    void write( void const* buf, size_t len ) { mFW.write(buf,len); }

    static const size_t BUF_SIZ = 16384;

    FileWriter mFW;
    char mBuf[BUF_SIZ];
    char* mpBuf;
};

/// Reader of binary streams.
class BinaryReader
{
public:
    BinaryReader( char const* filename, bool checkHeader = true )
    : mFR(filename), mpBuf(mBuf), mpEnd(mBuf)
    { if ( checkHeader ) testToken(); }

    /// Construct from something with a c_str() member (like string or String)
    template <class C>
    explicit BinaryReader( C const& filename, bool checkHeader=true,
                                char const*(C::*)() const=&C::c_str )
    : mFR(filename.c_str()), mpBuf(mBuf), mpEnd(mBuf)
    { if ( checkHeader ) testToken(); }

    BinaryReader( int fd, char const* pseudoFilename )
    : mFR(fd,pseudoFilename), mpBuf(mBuf), mpEnd(mBuf)
    {}

    // default destructor is OK

    std::string const& getFilename() const { return mFR.getFilename(); }

    /// Are we at the end of the stream?
    bool atEOF() { return mpBuf == mpEnd && !fillBuf(); }

    /// Read something.
    template <class T>
    T& read( T* pT )
    { readVal(pT,Serializability<T>()); return *pT; }

    /// Read an array.
    template <class T>
    void read( T* begin, T* end )
    { readArray(begin,end,Serializability<T>()); }

    /// Read a sequence delimited by an iterator pair.
    template <class Itr>
    void readItr( Itr begin, Itr const& end )
    { while ( begin != end ) { read(&*begin); ++begin; } }

    /// Read a single object from a binary file in one step.
    template <class T>
    static T& readFile( char const* filename, T* pObj )
    { BinaryReader rdr(filename);
      rdr.read(pObj);
      ForceAssert(rdr.atEOF());
      return *pObj; }

    /// Read a single object from a binary file in one step.
    /// The file is described by something with a c_str() member.
    template <class C, class T>
    static void readFile( C const& filename, T* obj,
                                char const*(C::*)() const=&C::c_str )
    { readFile(filename.c_str(),obj); }

    void seek( size_t pos )
    { mFR.seek(pos); mpBuf = mpEnd = mBuf; }

    template <class T>
    static size_t externalSizeof( T* arg )
    { return externalSizeof(arg,Serializability<T>()); }

private:
    BinaryReader( BinaryReader const& ); // unimplemented -- no copying
    BinaryReader& operator=( BinaryReader const& ); // unimplemented -- no copying

    template <class T>
    void readVal( T* pT, TriviallySerializable )
    { read(pT,sizeof(T)); }

    /// Read something that handles its own serialization.
    template <class T>
    void readVal( T* pVal, SelfSerializable )
    { pVal->readBinary(*this); }

    template <class T>
    void readVal( T* pT, ExternallySerializable )
    { readBinary(pT,*this); }

    template <class T>
    void readArray( T* begin, T* end, TriviallySerializable )
    { read(begin,sizeof(T)*(end-begin)); }

    template <class T>
    void readArray( T* begin, T* end, SelfSerializable )
    { while ( begin != end ) read(begin++); }

    template <class T>
    void readArray( T* begin, T* end, ExternallySerializable )
    { while ( begin != end ) read(begin++); }

    void read( void* pVal, size_t len )
    { size_t remain = mpEnd - mpBuf;
      if ( remain >= len )
      { memcpy(pVal,mpBuf,len); mpBuf += len; }
      else
      { if ( remain )
        { memcpy(pVal,mpBuf,remain);
          mpBuf = mpEnd;
          len -= remain;
          pVal = static_cast<char*>(pVal) + remain; }
        if ( len >= sizeof(mBuf) )
            mFR.read(pVal,len);
        else
            readLoop(static_cast<char*>(pVal),len); } }

    void readLoop( char* buf, size_t len );

    size_t fillBuf()
    { size_t result = mFR.readOnce(mBuf,sizeof(mBuf));
      mpBuf = mBuf; mpEnd = mBuf + result;
      return result; }

    template <class T>
    static size_t externalSizeof( T*, TriviallySerializable )
    { return sizeof(T); }

    template <class T>
    static size_t externalSizeof( T*, SelfSerializable )
    { return T::externalSizeof(); }

    template <class T>
    static size_t externalSizeof( T* arg, ExternallySerializable )
    { return serializedSizeof(arg); }

    void testToken();

    static size_t const BUF_SIZ = 16384;

    FileReader mFR;
    char mBuf[BUF_SIZ];
    char* mpBuf;
    char* mpEnd;
};

/// This is a very fragile class.  You parameterize it on some vector-like type,
/// and it allows you to read the elements of that vector from a binary stream
/// one by one without reconstituting the whole vector.
/// However:  It assumes that the vector was written to the stream as a
/// V::size_type followed by the V::value_type elements (as, indeed, all vectors
/// do, as of this writing).  If this isn't true, the whole thing will
/// likely explode in your face.
template <class V>
class BinaryIteratingReader
{
public:
    typedef typename V::size_type size_type;
    typedef typename V::value_type value_type;

    BinaryIteratingReader( char const* filename )
    : mRdr(filename)
    { mRdr.read(&mCount); }

    BinaryIteratingReader( std::istream& is )
    : mRdr(is)
    { mRdr.read(&mCount); }

    // compiler-supplied destructor is OK

    /// number of elements left to read
    size_type remaining() const { return mCount; }

    /// get the next element
    value_type next()
    { AssertGt(mCount,0u);
      --mCount;
      value_type result;
      mRdr.read(&result);
      return result; }

    class iterator
    : public std::iterator<std::input_iterator_tag,value_type>
    {
    public:
        iterator( BinaryIteratingReader* pIRdr = 0 )
        : mpIRdr(pIRdr)
        { operator++(); }

        // compiler-supplied destructor and copying are OK

        value_type const& operator*() const { Assert(mOK); return mVal; }
        value_type const* operator->() const { Assert(mOK); return &mVal; }

        iterator& operator++()
        { if ( (mOK = mpIRdr && mpIRdr->remaining()) ) mVal = mpIRdr->next();
          return *this; }

        iterator operator++( int )
        { iterator tmp(*this); operator++(); return tmp; }

        friend bool operator==( iterator const& itr1, iterator const& itr2 )
        { return itr1.mOK == itr2.mOK; }

        friend bool operator!=( iterator const& itr1, iterator const& itr2 )
        { return itr1.mOK != itr2.mOK; }

    private:
        BinaryIteratingReader* mpIRdr;
        value_type mVal;
        bool mOK; // indicates that mVal has a valid value
    };

    iterator begin() { return iterator(this); }
    iterator end() { return iterator(0); }

private:
    BinaryIteratingReader( BinaryIteratingReader const& ); // unimplemented -- no copying
    BinaryIteratingReader& operator=( BinaryIteratingReader const& ); // unimplemented -- no copying

    BinaryReader mRdr;
    typename V::size_type mCount;
};

/// Same, fragile idea as above.  Allow a vector-like type to be written
/// incrementally.  Closing patches the vector size, which is initially
/// recorded as 0.
template <class V>
class BinaryIterativeWriter : private BinaryWriter
{
public:
    BinaryIterativeWriter( char const* filename )
    : BinaryWriter(filename), mCount(0)
    { BinaryWriter::write(mCount); }

    BinaryIterativeWriter& write( typename V::value_type const& val )
    { BinaryWriter::write(val); mCount += 1; return *this; }

    void close()
    { BinaryWriter::patchSize(mCount); BinaryWriter::close(); }

private:
    typename V::size_type mCount;
};

/*
 * specializations for some commonly-used STL classes
 */
template <class T1, class T2>
struct Serializability<std::pair<T1,T2> > : public ExternallySerializable {};

/// Write a std::pair.
template <typename T1, typename T2>
inline size_t writeBinary( BinaryWriter& writer, std::pair<T1,T2> const& val )
{ size_t len = writer.write(val.first);
  return len+writer.write(val.second); }

/// Read a std::pair.
template <typename T1, typename T2>
inline void readBinary( std::pair<T1,T2>* pVal, BinaryReader& reader )
{ reader.read(&pVal->first); reader.read(&pVal->second); }

template <typename T1, typename T2>
inline size_t serializedSizeof( std::pair<T1,T2>* )
{ size_t sz1 = BinaryReader::externalSizeof(static_cast<T1*>(0));
  size_t sz2 = BinaryReader::externalSizeof(static_cast<T2*>(0));
  return sz1 && sz2 ? sz1+sz2 : 0UL; }

EXTERNALLY_SERIALIZABLE(std::string);

/// Write a std::string.
inline size_t writeBinary( BinaryWriter& writer, std::string const& str )
{ std::string::size_type nnn = str.size();
  size_t len = writer.write(nnn);
  if ( len )
  { char const* buf = &str[0];
    len += writer.write(buf,buf+nnn); }
  return len; }

/// Read a std::string.
inline void readBinary( std::string* pStr, BinaryReader& reader )
{ std::string::size_type sz; reader.read(&sz); pStr->resize(sz);
  if ( sz ) { char* buf = &(*pStr)[0]; reader.read(buf,buf+sz); } }

inline size_t serializedSizeof( std::string* )
{ return 0UL; }

template <class T, class A>
struct Serializability<std::vector<T,A> > : public ExternallySerializable {};

/// Write a std::vector
template <class T>
inline size_t writeBinary( BinaryWriter& writer, std::vector<T> const& vec )
{ typename std::vector<T>::size_type nnn = vec.size();
  size_t len = writer.write(nnn);
  T const* buf = &vec[0];
  return len+writer.write(buf,buf+nnn); }

/// Read a std::vector
template <class T>
inline void readBinary( std::vector<T>* pVec, BinaryReader& reader )
{ typename std::vector<T>::size_type sz; reader.read(&sz); pVec->resize(sz);
  T* buf = &(*pVec)[0]; reader.read(buf,buf+sz); }

template <class T>
inline size_t serializedSizeof( std::vector<T>* )
{ return 0UL; }

#endif /* FEUDAL_BINARYSTREAM_H_ */
