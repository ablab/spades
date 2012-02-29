///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/**
   Header: SemanticTypes.h

   Code for defining semantic types: different data types that may have the same
   physical representation but represent conceptually different and incomparable
   entities.  Semantic types allow such distinctions to be formally represented
   in the code.

   @file
*/

#ifndef __INCLUDE_SemanticTypes_h
#define __INCLUDE_SemanticTypes_h

/*
   Conditional define: SEMANTIC_TYPES_STRICT

   When defined, stricter compile-time checking is done for semantic types,
   at the possible expense of runtime efficiency (since we represent each
   semantic type as a separate class rather than as a built-in type for which
   the compiler and/or standard libraries may have specialized optimizations).
   When not defined, all <SemanticType()> definitions reduce to simple typedefs.

   So, the recommended usage is to define SEMANTIC_TYPES_STRICT during development,
   so you can catch more bugs, but to undefine it when building the production
   executable.
*/
#ifndef SEMANTIC_TYPES_STRICT

// Macro: SemanticType
// Defines a new semantic type, based on the given physical type.
// Right now it's just a typedef, but could be made to declare
// semanticType as a subclass of physicalType.
// For semantic types defined on built-in physical types (int, char etc) please use
// <SemanticTypeStd()> instead of this macro.
#define SemanticType(physicalType, semanticType) typedef physicalType semanticType

// Macro: SemanticTypeStd
// Defines a new semantic type, based on the given built-in physical type.
// Right now it's just a typedef, but could be made to declare
// semanticType as a subclass of physicalType.  Use this instead of
// <SemanticType()>.
#define SemanticTypeStd(physicalType, semanticType) typedef physicalType semanticType

#define IdxBy(t)

#else
// if SEMANTIC_TYPES_STRICT is defined

#error "Not implemented yet"

template <class T>
class SemType {
 public:
  typedef T value_type;
  
  T val;

  SemType() { }
  SemType(const T& _val): val(_val) { }

  SemType& operator=(const T& _val) { val = _val; return *this; }

  operator T() const { return val; }
};  // class SemType

#define SemanticTypeStd(physicalType, semanticType) class semanticType: public SemType<physicalType> {  \
     public:                                                                                            \
      typedef SemType<physicalType> PARENT;                                                             \
        semanticType() { }                                                                              \
        semanticType(const physicalType& _val): PARENT(_val) { }                                        \
        semanticType& operator=(const physicalType& _val) {                                             \
	  PARENT::operator=(_val); return *this;                                                        \
	}                                                                                               \
  }       

#define SemanticType(physicalType, semanticType) class semanticType: public physicalType {       \
       public:                                                                                   \
          typedef physicalType PARENT;                                                           \
          semanticType() { }                                                                     \
          semanticType( const physicalType& _val ): PARENT(_val) { }                             \
          semanticType& operator=( const physicalType& _val ) {                                  \
	    PARENT::operator=(_val); return *this;                                               \
	  }                                                                                      \
          operator PARENT() const { return PARENT(val); }                                        \
    }


#endif
// if SEMANTIC_TYPES_STRICT defined

#endif
// #ifndef __INCLUDE_SemanticTypes_h




