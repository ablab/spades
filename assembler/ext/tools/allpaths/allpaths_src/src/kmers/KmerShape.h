/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////


#ifndef __INCLUDE_KmerShape_h
#define __INCLUDE_KmerShape_h

#include "Basevector.h"
#include "CoreTools.h"
#include "ParseSet.h"
#include "CommonSemanticTypes.h"
#include "kmers/SupportedKmerShapes.h"

/**
   Type Concept: KmerShape

   Defines a shape for extracting kmers from <base vectors>: when
   extracting a kmer from a base vector starting at a given position,
   which positions in the base vector will be used?  By default, K
   consequtive positions are used, but other kmer shapes are possible.
   A kmer shape is defined by a list of offsets from the start of the
   kmer.  A particular kmer shape defines kmers of a particular size.
   There are several groups of kmer shapes -- the default kmer shape
   consisting of K consequtive bases, a kmer shape with a hole in the
   middle, a "zebra" kmer shape that extracts every other base; these
   we call <kmer forms>.  Typically, for each kmer form there is a
   templatized C++ class, whose instantiations define kmer shapes with
   that kmer form.  Each kmer shape is represented by a specific
   concrete C++ class which models the KmerShape concept.  (Note that
   there is no actual class named KmerShape; the documentation of
   class members below describes the members that a class must have in
   order to model the KmerShape concept.)  All kmer forms are
   parameterized on kmer size, and may have additional parameters.  A
   kmer form together with concrete values for all its parameters
   gives a particular kmer shape.

   To underline this again: the word "kmer shape" means a specific set
   of offsets, which always extracts kmers of a particular kmer size,
   rather than a general "shape" of kmers; the latter is referred to
   as kmer form.  Thus, a kmer shape is a particular instantiation of
   a kmer form.  A kmer form is a template whose instantiations give
   kmer shapes. Simply "default kmer form" is not a kmer shape;
   "default kmer form for K=8" is.  Kmer size is thus part of kmer
   shape; it is given by the <KSIZE> constant of the type representing
   the kmer shape.  When passing template arguments, it is therefore
   enough to pass the type representing the kmer shape; the kmer size
   is then available through the KSIZE constant.

   Classes that model the KmerShape concept and represent kmer shapes,
   define a method <extractKmer()> which extracts the kmer of the
   given shape from a location on a base vector.

   Programmatically, kmer shapes are identified by a <KmerShapeId>,
   which can be constructed from a string.  See <KmerShapeId(const
   String&)> for documentation on which strings denote which kmer
   shapes.  The string "K" where K is kmer size will always denote the
   default kmer shape of K consequtive bases.  The set of kmer shapes
   supported by the code is controlled by definitions in
   <SupportedKmerShapes.h>.

   To invoke the right instantiation of code templatized on kmer shape, you do the following:

        - take a string command-line argument identifying the shape,
	
	- define a macro that takes kmer shape (i.e. the name of the
        concrete class that models KmerShape and denotes a particular
        kmer shape) and calls your code instantiated with that kmer
        shape,
	
	- invoke the <DISPATCH_ON_KSHAPE()> macro, passing it the
          command-line argument and your macro as arguments.
	
	- define a macro that takes two arguments, kmer shape and a
	  second value (normally ignored), and explicitly instantiates
	  your templatized code for the passed kmer shape
	  
	- in the declarations section of your .cc file, invoke the
	<FOR_ALL_KSHAPE()> macro passing it your instantiation macro
	as first argument, and a dummy value as second argument.  (The
	value you pass as second argument to <FOR_ALL_KSHAPE()> will
	be passed as second argument to each invocation of your
	explicit instantiation macro.)

   See <FindStrongKmers()> and <WriteKmerFrequencies()> for an example.

   If your code is parameterized by kmer size only, but not by kmer
   shape (for example it works with kmers that have already been
   extracted, and does not care about their original shape), then you
   can use the <DISPATCH_ON_K()> and <FOR_ALL_K()> macros instead of
   <DISPATCH_ON_KSHAPE()> and <FOR_ALL_KSHAPE()>, respectively.

   Some of the supported kmer forms:

        the default kmer form - just the normal kmer shape, extracting
	   K consequtive kmers.  Use <KmerShapeDefaultClass()> to
	   denote kmer shapes of this form.  On the command line, use
	   the string "K" (as in "16") to specify kmers of size K of
	   the default form.
	   
	the midgap kmer form - has a gap of G bases in the middle.
	   Extracts K/2 bases followed by a G-base gap followed by
	   extracting the remaining K/2 bases.  Use
	   <KmerShapeMidGapClass()> to denote kmer shapes of this
	   form.  On the command line, use the string "KgG" (as in
	   "16g4") to specify kmers of size K with a gap of G bases in
	   the middle.
	   
	the zebra shape - extracts every other base, until K bases are
	   extrated.  Use <KmerShapeZebraClass()> to denote kmer
	   shapes of this form.  On the command line, use the string
	   "Kz" to denote kmers of zebra form of size K.
   
   Note that for a kmer, exactly K bases are always extracted from the
   base vector; and after extraction, the k bases are represented as a
   continuous k-base array.  However, when we have a <kmer_record>
   representing the occurrence of a kmer on a read, exactly what the
   occurrence means depends on the kmer shape.

   *NOTE*: the shape _must_ be symmetric -- that is, if flipped about
   the center, we must get the same shape!  This is required so that
   the reverse complement of the extracted kmer corresponds to the
   kmer extracted by the same kmer shape from the complementary
   strand.

   Note that we also assume that the kmer size K is even.  This is
   normally true, as <SORT_KMERS()> requires it to be divisible by
   four.

   No instances of a KmerShape class are ever created, but the static
   members of a KmerShape class are used to implement a particular
   version of kmer extraction.

  A class that is a model of KmerShape must define the _static_ members described below. 

  Constant: KSIZE
  The size of the kmer: the number of bases extracted from the base vector.
  >static const int KSIZE;

  Constant: KSPAN
  The span from which the kmer is extracted: the distance between the largest and smallest offset
  from the start of the kmer extraction region.
  >static const int KSPAN;

  Method: getKmerSize
  Return the kmer size
  >static int getKmerSize();

  Method: getShapeOffset
  Return the offset of the i'th base of the shape, from the leftmost base of the shape.
  >static unsigned int getShapeOffset(int posInShape);

  Method: getNumGaps
  Return the number of inner gaps (where each position counts as a separate gap) in the shape.
  >static unsigned int getNumGaps();

  Method: getGapOffset
  Return the offset of the i'th gap, from the leftmost base of the shape.
  >static unsigned int getGapOffset(int gapNum);

  Method: getSpan
  The size of the region from which we gather a shape.
  >static unsigned int getSpan() { return K + GAPLEN; }

  Method: extractKmer
  
  Extract the kmer from the given position of a read, and store it in a kmer.
  The kmer shape *must* be symmetric!

  Parameters:
  
      fwdKmer -the extracted kmer (forward version) is put here
      read - the read from which to extract the kmer
      posInRead - the position in the read, at which the kmer begins.

  See also <SORT_CORE>, <SupportedKmerShapes.h>

  Of course, the methods <getShapeOffset()> and <getGapOffset()>
  already give enough information to implement kmer extraction;
  nevertheless, each kmer shape should explicitly define its own
  extractKmer() method specific to the kmer form, which is often more
  efficient than a generic method that worked for all kmer forms would
  be.

  >static void extractKmer( basevector& fwdKmer, const basevector& read, int posInRead );

  Method: getId

  Get the id denoting this kmer shape.   See <KmerShapeId>.
*/

/// Local function: KmerShapeDefaultStringId
/// Return the string id denoting kmers of the <default form> and given size --
/// must be just the string version of the size, for backwards compatibility.
inline String KmerShapeDefaultStringId(int K) { return ToString(K); }

/**
   Class: KmerShapeId

   A value that identifies/represents a particular <kmer shape>.  This
   is for when we need to pass around kmer shapes as a normal function
   argument and treat it as a normal value, rather than passing it
   around as a type in a template argument.

   Returned by <KmerShape::getId()>.

   See also: <DISPATCH_ON_KSHAPE()>.
*/
class KmerShapeId {
 public:
  /// Constructor: KmerShapeId(int)
  /// Create a kmer shape id for the <default kmer form> of the given size.
  explicit KmerShapeId(int ksize): id_(KmerShapeDefaultStringId(ksize)), ksize_(ksize) {
  }
  
  /// Constructor: KmerShapeId(const String&)
  /// Create a kmer shape id from a string description.
  //  The currently supported kmer shape ids are:
  //
  //        K - the default kmer form of size K
  //        KgG - kmer of size K with gap of size G in the middle
  //        Kz - "zebra" kmer of size K, extracting every other base.
  explicit KmerShapeId(const String& id): id_(id) {
    // check that this id is a valid designator.
    extractKmerSize_();
  }

  KmerShapeId(const KmerShapeId& src):
    id_(src.id_), ksize_(src.ksize_) { }

  KmerShapeId(): id_(""), ksize_(0) { }

  KmerShapeId& operator= (const KmerShapeId& ksi) {
    id_ = ksi.id_; ksize_ = ksi.ksize_; return *this;
  }

  KmerShapeId& operator= (const String& ksi) {
    id_ = ksi;  extractKmerSize_(); return *this;
  }

  KmerShapeId& operator= (int ksize) {
    id_ = KmerShapeDefaultStringId(ksize); ksize_ = ksize; return *this;
  }

  //operator int() const { return ksize_; }

 private:
  /// Field: id_
  /// A unique string identifying this kmer shape, as returned by <KmerShape::getStringId()>.
  String id_;

  /// Field: ksize_
  /// The kmer size of the kmer shape represented here.
  int ksize_;

  /// Private method: extractKmerSize_
  /// Determine <ksize_> from <id_>.
  void extractKmerSize_() {
    // for now, the convention is that the string starts with an integer identifying the kmer size.
    String num;
    for (int i = 0; i < id_.isize() && isdigit(id_[i]); i++)
      num += id_[i];

    ksize_ = num.Int();
  }
  
  friend bool operator==(const KmerShapeId& ksi1, const KmerShapeId& ksi2);
  friend bool operator<(const KmerShapeId& ksi1, const KmerShapeId& ksi2);
  
  friend int GetKmerSize(const KmerShapeId& kmerShapeStringId);
  friend String ToString(const KmerShapeId& ksi);
  friend ostream& operator<< ( ostream& out, const KmerShapeId& ksi );
};  // class KmerShapeId

inline bool operator==(const KmerShapeId& ksi1, const KmerShapeId& ksi2) {
  return ksi1.id_ == ksi2.id_;
}

inline bool operator<(const KmerShapeId& ksi1, const KmerShapeId& ksi2) {
  return ksi1.id_ < ksi2.id_;
}

inline bool operator!=(const KmerShapeId& ksi1, const KmerShapeId& ksi2) {
  return !(ksi1 == ksi2);
}


inline bool operator>(const KmerShapeId& ksi1, const KmerShapeId& ksi2) {
  return !( ksi1 < ksi2  ||  ksi1 == ksi2 );
}

inline bool operator>=(const KmerShapeId& ksi1, const KmerShapeId& ksi2) {
  return !( ksi1 < ksi2 );
}

inline bool operator<=(const KmerShapeId& ksi1, const KmerShapeId& ksi2) {
  return !( ksi1 > ksi2 );
}


inline ostream& operator<< ( ostream& out, const KmerShapeId& ksi ) {
  out << ToString(ksi);
  return out;
}

/// Function: GetKmerSize
/// Return the size of the kmers represented by the kmer shape identified
/// by the given <KmerShapeId>.
inline int GetKmerSize(const KmerShapeId& kmerShapeId) {
  return kmerShapeId.ksize_;
}

/**
   Function: ParseKmerShapeIdSet

   Given a string describing a set of <kmer shape ids>, parse it into a
   vec of <KmerShapeId>.  See <ParseStringSet()> for explanation of the
   string format.
*/
inline void ParseKmerShapeIdSet( String descrip, vec<KmerShapeId>& answer,
				 Bool ABORT_IF_BAD = False, bool sortAnswer = true  ) {
  vec<String> answerStrings;
  ParseStringSet( descrip, answerStrings );
  if (sortAnswer) {
    Sort( answerStrings );
  }
  for (int i=0; i < answerStrings.isize(); i++)
    answer.push_back( KmerShapeId( answerStrings[i] ) );
}

inline String ToString(const KmerShapeId& ksi) { return ksi.id_; }

#define CommandArgument_KShape(NAME) \
     KmerShapeId NAME; \
     if ( command.GetHelpOnly() ) \
       command.AddArgHelp( #NAME, "kshape", "", "<required>" ); \
     else \
       NAME = command.GetStringValue( #NAME )

#define CommandArgument_KShape_OrDefault(NAME, DEFAULT) \
     KmerShapeId NAME; \
     if ( command.GetHelpOnly() ) \
       command.AddArgHelp( #NAME, "kshape", "", DEFAULT ); \
     else \
       NAME = command.GetStringValue( #NAME, "", DEFAULT );

#define CommandArgument_KShapes(NAME) \
     vec<KmerShapeId> NAME; \
     if ( command.GetHelpOnly() ) \
       command.AddArgHelp( #NAME, "kshapes", "", "<required>" ); \
     else \
       ParseKmerShapeIdSet( command.GetStringValue( #NAME ), NAME )

#define CommandArgument_KShapes_OrDefault(NAME, DEFAULT) \
     vec<KmerShapeId> NAME; \
     if ( command.GetHelpOnly() ) \
       command.AddArgHelp( #NAME, "kshapes", "", DEFAULT ); \
     else \
       ParseKmerShapeIdSet(command.GetStringValue( #NAME, "", DEFAULT ), NAME )


#define CommandArgument_KShapes2(NAME, NAME2) \
     vec<KmerShapeId> NAME2; \
     if ( command.GetHelpOnly() ) \
       command.AddArgHelp( #NAME, "kshapes", "", "<required>" ); \
     else \
       ParseKmerShapeIdSet( command.GetStringValue( #NAME ), NAME2 )


/**
   Class: kmer_shape_mid_gap
   
   The kmer shape with one specified gap of the given length in the middle.
   Make sure that the gap is not so large as to make the <span> exceed the typical read length!!
*/
template <int K, int GAPLEN=0>
class kmer_shape_mid_gap {
 public:
  /// Constant: KSIZE
  /// The size of the kmer.
  static const int KSIZE = K;

  /// Constant: KSPAN
  /// The span from which the kmer is extracted
  static const int KSPAN = K + GAPLEN;
  
  /// Return the kmer size.
  static int getKmerSize() { return K; }
  
  /// Return the offset of the i'th base of the shape, from the leftmost base of the shape.
  static unsigned int getShapeOffset(int posInShape) { return posInShape < K/2 ? posInShape : posInShape+GAPLEN; }

  /// Return the number of inner gaps (where each position counts as a separate gap) in the shape.
  static unsigned int getNumGaps() { return GAPLEN; }

  /// Return the offset of the i'th gap, from the leftmost base of the shape.
  static unsigned int getGapOffset(int gapNum) {
    Assert( 0 < gapNum  &&  gapNum < GAPLEN );
    return K/2 + gapNum;
  }
  
  /// The size of the region from which we gather a shape.
  static unsigned int getSpan() { return K + GAPLEN; }

  /**
     Method: extractKmer
     
     Extract the kmer from the given position of a read, and store it in a kmer.
     
     Input parameters:
     
        read - the read from which to extract the kmer
        posInRead - the position in the read, at which the kmer begins.

     Output parameters:
     
        fwdKmer - the extracted kmer (forward version) is put here

     See also: <SORT_CORE()>.
  */
  static void extractKmer( basevector& fwdKmer, const basevector& read, int posInRead ) {
    fwdKmer.SetToSubOf(read, posInRead, K);
    CopyBases(read, posInRead+K/2+GAPLEN, fwdKmer,
	      K/2 /* start writing in the middle of extractedKmer*/,
	      K/2 /* copy the remainder */);
  }

 // Method: getId
 // Return a <KmerShapeId> object that uniquely identifies this particular kmer shape
  static KmerShapeId getId() { return KmerShapeId( getStringId() ); }

 private:
 //  Method: getStringId
 //  Return a string that uniquely identifies this particular kmer shape.
 static String getStringId() { return ToString(K) + "g" + ToString(GAPLEN); }
 
};  // class kmer_shape_mid_gap 

/**
   The default kmer shape, ungapped: K contiguous letters.
   Defined as a specialization of kmer_shape_mid_gap<K,0> because this allows for uniform processing in some cases
   (for example, see FOR_ALL_K_GAP() in MacroUtils.h).
*/
template <int K>
class kmer_shape_mid_gap<K,0> {
public:
  // Constant: KSIZE
  // The size of the kmer.
  static const int KSIZE = K;

  // Constant: KSPAN
  // The span from which the kmer is extracted
  static const int KSPAN = K;
  
  
  /// Return the kmer size.
  static int getKmerSize() { return K; }
  
  /// Return the offset of the i'th base of the shape.
  static unsigned int getShapeOffset(int posInShape) { return posInShape; }
  
  /// Return the number of inner gaps (where each position counts as a separate gap) in the shape.
  static unsigned int getNumGaps() { return 0; }
  
  /// Return the offset of the i'th gap, from the leftmost base of the shape.
  static unsigned int getGapOffset(int gapNum) { ForceAssert( 0 ) ; return 0; }
    
  /// The size of the region from which we gather a shape.
  static unsigned int getSpan() { return K; }

  /**
     Method: extractKmer
     
     Extract the shape from the given position of a read, and store it in a kmer.

     Parameters:

      fwdKmer - the extracted kmer (forward version) is put here
      read - the read from which to extract the kmer
      posInRead - the position in the read, at which the kmer begins.

     See also: <SORT_CORE>
  */
  static void extractKmer( basevector& fwdKmer, const basevector& read, int posInRead ) {
    fwdKmer.SetToSubOf(read, posInRead, K);
  }

 // Method: getId
 // Return a <KmerShapeId> object that uniquely identifies this particular kmer shape
 static KmerShapeId getId() { return KmerShapeId( getStringId() ); }

 private:
  // Method: getStringId
  // Return a string that uniquely identifies this particular kmer shape.
  static String getStringId() { return KmerShapeDefaultStringId(K); }
 
};  // class kmer_shape_mid_gap<K,0>

/**
   Class: kmer_shape_zebra

   Create a "zebra" kmer shape where shape positions and gap positions alternate.
   Like this:
   > XOXOXOXOOXOXOXOX
   In the above, the X's represent shape positions while the O's represent gaps in the shape.
   Remember that the shape must start and end with a shape position rather than a gap,
   that the shape must be symmetric, and the kmer size (number of shape positions) must be
   divisible by four (this last requirement is imposed in particular by the current (04-12-07)
   implementation of <SortKmers()>).
*/
template <int K>
class kmer_shape_zebra {
public:
  /// Constant: KSIZE
  /// The size of the kmer.
  static const int KSIZE = K;

  /// Constant: KSPAN
  /// The span from which the kmer is extracted
  static const int KSPAN = 2 * K;
  
  /// Return the kmer size.
  static int getKmerSize() { return K; }
  
  /// Return the offset of the i'th base of the shape.
  static unsigned int getShapeOffset(int posInShape) { return 2*posInShape + (posInShape < K/2  ?  0 : 1) ; }
  
  /// Return the number of inner gaps (where each position counts as a separate gap) in the shape.
  static unsigned int getNumGaps() { return K; }
  
  /// Return the offset of the i'th gap, from the leftmost base of the shape.
  static unsigned int getGapOffset(int gapNum) { return 2*gapNum + (gapNum < K/2 ? 1 : 0); }
    
  /// The size of the region from which we gather a shape.
  static unsigned int getSpan() { return KSPAN; }

  /**
     Extract the shape from the given position of a read, and store it in a kmer.

     Parameters:

       fwdKmer - the extracted kmer (forward version) is put here
       read - the read from which to extract the kmer
       posInRead - the position in the read, at which the kmer begins.

     See also: <SORT_CORE>.

  */
  static void extractKmer( basevector& fwdKmer, const basevector& read, int posInRead ) {
    // Inefficiency: Naive extraction of zebra kmers
    // It might be possible to optimize this method by using a precomputed lookup table.
    int i, readPos;
    for (i = 0, readPos = 0; i < K/2; i++, readPos += 2)
      fwdKmer.Set( i, read[ readPos ] );
    for (readPos++; i < K; i++, readPos += 2)
      fwdKmer.Set( i, read[ readPos ] );
  }

 // Method: getId
 // Return a <KmerShapeId> object that uniquely identifies this particular kmer shape
 static KmerShapeId getId() { return getStringId(); }

 private:
  // Method: getStringId
  // Return a string that uniquely identifies this particular kmer shape.
  static String getStringId() { return ToString(K) + "z"; }
 
};  // class kmer_shape_zebra<K>


/**
   Class: KmerShapeDefault

   The standard kmer shape, which extracts K consequtive bases: a kmer shape without any gaps.
  
   Alternative name for kmer_shape_mid_gap<K,0>, designating the default kmer shape.

   The definition below works because the second template argument of kmer_shape_mid_gap
   defaults to zero, and the template specialization
   kmer_shape_mid_gap<K,0> of kmer_shape_mid_gap<K,GAP> implements the default kmer shape.
*/
#define KmerShapeDefault kmer_shape_mid_gap

#define KmerShapeDefaultClass(K) KmerShapeDefault<K,0>

/*
   Macro: KmerShapeMidGapType
   
   Return the typedef name for kmer_shape_mid_gap<K,GAP>.  The resulting name has no commas, which makes it
   usable in various macros -- if you use kmer_shape_mid_gap<K,GAP> as an actual macro argument, the preprocessor
   will get confused, thinking it got two tokens separated by a comma, since it does not treat angle brackets as paired brackets.
*/
#define KmerShapeMidGapType(K,GAP) kmer_shape_mid_gap_ ## K ## _ ## GAP
#define KmerShapeDefaultType(K) KmerShapeMidGapType(K,0)

template < int K >
struct KmerShapeDflt {
  typedef KmerShapeDefault< K, 0 > type;
};

/**
   Create a typedef name for a particular instantiation of kmer_shape_mid_gap<K,GAP>.   THese are useful
   in situations where you need  a single identifier (with no commas) denoting the type -- for example,
   when you want to pass it in a macro argument.
 */
#define CreateKmerShapeMidGapTypedef(K,GAP) typedef kmer_shape_mid_gap<K,GAP> KmerShapeMidGapType(K,GAP)

FOR_ALL_K(CreateKmerShapeMidGapTypedef, 0);
FOR_ALL_K(CreateKmerShapeMidGapTypedef, 1);
FOR_ALL_K(CreateKmerShapeMidGapTypedef, 4);
FOR_ALL_K(CreateKmerShapeMidGapTypedef, 8);

/// Macro: KmerShapeZebraClass
/// Return the class denoting kmers of the <zebra form> and given kmer size.
#define KmerShapeZebraClass(K) kmer_shape_zebra<K>


/// Local macro: KSHAPE_CASE_GET_shapeId
/// Given a (shape id, shape id handler macro) pair, extract the shape id.
/// Used by <KSHAPE_CASE()> which is used by <DISPATCH_ON_KSHAPE()>.
#define KSHAPE_CASE_GET_shapeId(shapeId, handleShape) shapeId

/// Local macro: KSHAPE_CASE_GET_handleShape
/// Given a (shape id, shape id handler macro) pair, extract the shape id handler
/// macro.
/// Used by <KSHAPE_CASE()> which is used by <DISPATCH_ON_KSHAPE()>.
#define KSHAPE_CASE_GET_handleShape(shapeId, handleShape) handleShape


/// Local macro: KSHAPE_CASE
/// Used by <DISPATCH_ON_KSHAPE()> to implement one case of the dispatcher,
/// handling a particular kmer shape.
#define KSHAPE_CASE(shapeType,shapeId_handleShape)                           \
   if ( KmerShapeId( KSHAPE_CASE_GET_shapeId shapeId_handleShape ) == shapeType::getId()) {  \
      foundKShape = true;                                                    \
      KSHAPE_CASE_GET_handleShape shapeId_handleShape (shapeType) ;          \
    }                                               
   
/**
   Macro: DISPATCH_ON_KSHAPE

   Call the given macro on the shape specified by the given <shape id>.
   
   Parameters:
   
      shapeId - the <shape id> value, of type KmerShapeId; can be a variable (not a constant)
      handleShape - a macro that takes one argument -- the <kmer shape> (the name of the class
         representing that shape), and specifies what to do for a given shape.
        
*/
#define DISPATCH_ON_KSHAPE(shapeId,handleShape) do {                        \
     bool foundKShape = false;                                              \
     FOR_ALL_KSHAPES(KSHAPE_CASE, (shapeId, handleShape));                  \
     if (!foundKShape) {                                                    \
         cout << "Not implemented for KSHAPE=" << shapeId << "." << endl;   \
         TracebackThisProcess();                                            \
     }                                                                      \
  } while(0)


#define DISPATCH_ON_MAIN_KSHAPE(shapeId,handleShape) do {                        \
     bool foundKShape = false;                                              \
     FOR_MAIN_KSHAPES(KSHAPE_CASE, (shapeId, handleShape));                  \
     if (!foundKShape) {                                                    \
         cout << "Not implemented for KSHAPE=" << shapeId << "." << endl;   \
         TracebackThisProcess();                                            \
     }                                                                      \
  } while(0)


/// Macro: NO_SHAPES
/// Use in <SupportedKmerShapes.h> when defining the <FOR_ALL_KSHAPES_A()> and similar
/// macros, if you want to leave a particular one of these macros empty.
/// This macro makes it ok to add a semicolon after instantiating these macros
/// in the declarations section of a .cc file.
#define NO_KSHAPES(X) typedef int no_kshapes_ ## X ## _ 

/*
   Macro: FOR_ALL_KSHAPES

   Call the specified macro for the <kmer shape set> that we want to support,
   passing the specified argument to the macro in addition to the
   kmer shape.  (If you need to pass more than one argument, see
   the definition of <DISPATCH_ON_KSHAPE()> for an illustration.)

   To define the list shapes for which this macro is called, define
   <FOR_ALL_KSHAPES_A()> and similar macros in <SupportedKmerShapes.h>.

   Use this macro to do explicit instantiation of code templatized on kmer shape.
   Because <DISPATCH_ON_K()> also uses this macro, when you use DISPATCH_ON_K()
   it will dispatch only on the instantiated kmer shapes, ensuring that the
   build goes through.

   See also <FOR_ALL_K()>.
*/
#define FOR_ALL_KSHAPES(M, arg) \
   FOR_ALL_KSHAPES_A(M, arg); \
   FOR_ALL_KSHAPES_B(M, arg); \
   FOR_ALL_KSHAPES_C(M, arg); \
   FOR_ALL_KSHAPES_D(M, arg)

#define FOR_MAIN_KSHAPES(M, arg) \
   M(KmerShapeDefaultType(20), arg); \
   M(KmerShapeDefaultType(21), arg) 

#define K_CASE_GET_KVAR(Kvar,handleK) Kvar
#define K_CASE_GET_handleK(Kvar,handleK) handleK

#define K_CASE(K, Kvar_handleK)  \
   if (K_CASE_GET_KVAR Kvar_handleK == K) {             \
      foundK = true;                                    \
      K_CASE_GET_handleK Kvar_handleK (K) ;             \
    }                                               
   
/**
   Macro: DISPATCH_ON_K_WITH_K_PLUS_1

   Take kmer size passed as a variable, and a macro that takes kmer size as a constant argument,
   and construct a dispatcher that calls the macro's body for the kmer size matching the value
   of the variable. Only allowed for values of K for which there is also a valid K+1 value.
   See FOR_ALL_K_WITH_K_PLUS_1
   
   Parameters:
   
      Kvar - the kmer size value; can be a variable (not a constant)
      handleK - a macro that takes one argument -- the kmer size, and specifies what to
         do for that kmer size.
        
*/
#define DISPATCH_ON_K_WITH_K_PLUS_1(Kvar,handleK) do {                               \
     bool foundK = false;                                              \
     FOR_ALL_K_WITH_K_PLUS_1(K_CASE, (Kvar, handleK));                               \
     if (!foundK) {                                                    \
         cout << "Not implemented for K=" << Kvar << "." << endl;      \
         TracebackThisProcess();                                       \
     }                                                                 \
  } while(0)


/**
   Macro: DISPATCH_ON_K

   Take kmer size passed as a variable, and a macro that takes kmer size as a constant argument,
   and construct a dispatcher that calls the macro's body for the kmer size matching the value
   of the variable.
   
   Parameters:
   
      Kvar - the kmer size value; can be a variable (not a constant)
      handleK - a macro that takes one argument -- the kmer size, and specifies what to
         do for that kmer size.
        
*/
#define DISPATCH_ON_K(Kvar,handleK) do {                               \
     bool foundK = false;                                              \
     FOR_ALL_K(K_CASE, (Kvar, handleK));                               \
     if (!foundK) {                                                    \
         cout << "Not implemented for K=" << Kvar << "." << endl;      \
         TracebackThisProcess();                                       \
     }                                                                 \
  } while(0)


/**
    Macro: FOR_SUPPORTED_K

    Take a constant kmer size value, and invoke code to handle that value.
*/
#define FOR_SUPPORTED_K(Kval,handleK) DISPATCH_ON_K(Kval, handleK)   

inline void ForceAssertSupportedK(int K) {
#define CHK_K(_K)  // do nothing: if K is not supported, we'll get an error message
  DISPATCH_ON_K(K, CHK_K);
}

inline void ForceAssertSupportedKShape(const KmerShapeId& ksi) {
#define CHK_KSHAPE(_KSHAPE)  // do nothing: if K is not supported, we'll get an error message
  DISPATCH_ON_KSHAPE(ksi, CHK_KSHAPE);
}

inline void ForceAssertSupportedKShapes(const vec<KmerShapeId>& ksis) {
  for (int i=0; i<ksis.isize(); i++)
    ForceAssertSupportedKShape(ksis[i]);
}

#endif
// #ifndef __INCLUDE_KmerShape_h

/*
   Term: kmer form

   A particular template for <kmer shapes>: for example, one form gives the default kmer shapes; another form
   gives kmer shapes with one hole in the middle; another form gives kmer shapes that extract every other base
   in a zebra-like fashion; and so on.

   All kmer forms are parameterized on kmer size, and
   may have additional parameters.  A kmer form together with concrete values for all its parameters gives a particular
   kmer shape.
*/

// Synonyms: Various synonyms
//   kmer shape - See <KmerShape>
//   shape of kmers - See <KmerShape>
//   shape - See <KmerShape>
//   comb - See <kmer shape>
//   KmerShape.h - See <KmerShape>



