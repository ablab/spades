///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file BinaryStreamTraits.h
 * \author tsharpe
 * \date Oct 16, 2009
 *
 * \brief
 */
#ifndef FEUDAL_BINARYSTREAMTRAITS_H_
#define FEUDAL_BINARYSTREAMTRAITS_H_

/// This is a trait that describes how classes should be serialized.
/// There are three categories:
///
/// TRIVIALLY_SERIALIZABLE
/// Easiest is trivially serializable:  classes that are trivially serializable
/// just read and write their bits, just as they are laid out in memory.  To be
/// eligible for trivial serialization your class must have only plain-old data
/// (POD) as instance data.  PODs are primitive types like int and double,
/// classes that are themselves trivially serializable, and fixed-length arrays
/// of all these. It must not contain references or pointers or classes (like
/// String or string, to choose a common example) that are not themselves
/// trivially serializable.  This is the easiest category to implement since you
/// need do nothing except to mark them as trivially serializable -- the reader
/// and writer handle their serialization for you.
/// Tell the system that your class is trivially serializable by using the macro
/// TRIVIALLY_SERIALIZABLE(myClass);
///
/// SELF_SERIALIZABLE
/// The next category is classes that take responsibility for their own
/// serialization.  A class that is self-serializable implements the following
/// three methods:
///    size_t writeBinary( BinaryWriter& writer ) const;
///    void readBinary( BinaryReader& reader );
///    static size_t externalSizeof();
/// The last method is only used in reading feudal files (not binary streams),
/// and is a hint that can make things more efficient.  If, in serializing your
/// class, you always write exactly the same number of bytes, then return that
/// number as your external size.  Otherwise return 0.
/// Tell the system that your class is self-serializable by using the macro
/// SELF_SERIALIZABLE(myClass);
///
/// EXTERNALLY_SERIALIZABLE
/// The last category is ExternallySerializable.  Externally serializable
/// classes rely on two helper functions:
/// size_t writeBinary( BinaryWriter&, T const& );
/// void readBinary( T*, BinaryReader& );
/// This is useful in the case where the class is some STL class, or something
/// else you can't easily alter.  (Self serialization is a better choice for
/// classes that you can alter.)  It provides for completely unintrusive
/// serialization, but is obviously only useful if your class has methods that
/// allow you to inspect its entire state from the "outside".
///
/// If you fail to declare a Serializability trait for some class which you then
/// try to serialize, you'll get a compile-time error like this:
/// error: invalid use of incomplete type 'struct Serializability<someClass>'
///
struct TriviallySerializable {};
struct SelfSerializable {};
struct ExternallySerializable {};

template <class T> struct Serializability;

#define TRIVIALLY_SERIALIZABLE(T) \
    template <> struct Serializability<T> : public TriviallySerializable {}
#define SELF_SERIALIZABLE(T) \
    template <> struct Serializability<T> : public SelfSerializable {}
#define EXTERNALLY_SERIALIZABLE(T) \
    template <> struct Serializability<T> : public ExternallySerializable {}

TRIVIALLY_SERIALIZABLE(bool);
TRIVIALLY_SERIALIZABLE(char);
TRIVIALLY_SERIALIZABLE(signed char);
TRIVIALLY_SERIALIZABLE(unsigned char);
TRIVIALLY_SERIALIZABLE(short);
TRIVIALLY_SERIALIZABLE(unsigned short);
TRIVIALLY_SERIALIZABLE(int);
TRIVIALLY_SERIALIZABLE(unsigned int);
TRIVIALLY_SERIALIZABLE(long);
TRIVIALLY_SERIALIZABLE(unsigned long);
TRIVIALLY_SERIALIZABLE(float);
TRIVIALLY_SERIALIZABLE(double);

#endif /* FEUDAL_BINARYSTREAMTRAITS_H_ */
