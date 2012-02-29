#ifndef __INCLUDE_SupportedKmerShapes_h
#define __INCLUDE_SupportedKmerShapes_h

/**
   Header file: SupportedKmerShapes.h

   Defines the set of <kmer shapes> supported by the system.  Change this header file
   to trim the set to just the kmer shapes you use for faster compilation, or
   to include a complete set of kmer shapes for production distribution.

   The supported set of kmer shapes is defined by a higher-order macro <FOR_ALL_KSHAPE()> that takes a
   macro as an argument and calls it for all supported kmer shapes.  This lets
   us explicitly instantiate all template classes and functions templatized
   on kmer shape, for the supported set of kmer shapes.  For each template
   that needs to be instantiated for all kmer shapes, you define a macro
   that instantiates the template for one kmer shape (passed as an argument),
   and call FOR_ALL_KSHAPES() passing it your macro as first argument.
   Your macro then gets invoked for each of the supported kmer shapes, and explicitly
   instantiates your template for each shape.

   To speed up compilation, FOR_ALL_KSHAPES() is actually broken up into several macros,
   FOR_ALL_KSHAPES_A through FOR_ALL_KSHAPES_I.  You can spread the supported kmer shapes
   among the macros as you like, but try to spread them evenly for best compilation time.

   For each kmer shape you want the system to support, you add the line
   >   M(KSHAPE, arg); \
   to the definition of one of the FOR_ALL_KSHAPES_[A-I] macros below,
   where KSHAPE is the shape you want to add.
   
   Please use the following macro invocations in place of KSHAPE.  These macros are defined
   in <KmerShape.h>.

      KmerShapeDefaultType(K) - the default kmer shape, K consequtive bases
      KmerShapeMidGapType(K,G) - the kmer shape with one G-base gap in the middle:
          K/2 bases followed by a G-long gap followed by K/2 bases.
      KmerShapeZebraType(K) - the kmer shape of alternating bases and gaps: XOXOOXOX for K=4,
         where the X's denote bases of the shape and O's denote gaps in the shape.

   To choose kmer shapes on the command line, use
      K - to specify the <default kmer shape> of size K
      KgG - to specify the <midgap kmer shape> of size K with gap of size G
      Kz - to specify the <zebra kmer shape> of size K

   In addition to adding the shapes you want the system to support to the definitions of the
   FOR_ALL_KSHAPES_[A-I] macros, you need add the kmer sizes of these kmer shapes to
   the <FOR_ALL_K()> macro.  Each size must be added exactly once, regardless of how many
   kmer shapes with that size are supported.
*/

#define DEFAULT_KMER_SHAPE_SET
#ifdef DEFAULT_KMER_SHAPE_SET

#define FOR_ALL_KSHAPES_A(M, arg) \
   M(KmerShapeDefaultType(4), arg); \
   M(KmerShapeDefaultType(8), arg); \
   M(KmerShapeDefaultType(12), arg); \
   M(KmerShapeDefaultType(16), arg); \
   M(KmerShapeDefaultType(19), arg); \
   M(KmerShapeDefaultType(20), arg); \
   M(KmerShapeDefaultType(21), arg); \
   M(KmerShapeDefaultType(24), arg); \
   M(KmerShapeDefaultType(26), arg)

#define FOR_ALL_KSHAPES_B(M, arg) \
   M(KmerShapeDefaultType(28), arg); \
   M(KmerShapeDefaultType(29), arg); \
   M(KmerShapeDefaultType(32), arg); \
   M(KmerShapeDefaultType(36), arg); \
   M(KmerShapeDefaultType(40), arg); \
   M(KmerShapeDefaultType(48), arg); \
   M(KmerShapeDefaultType(49), arg); \
   M(KmerShapeDefaultType(52), arg)

#define FOR_ALL_KSHAPES_C(M, arg) \
   M(KmerShapeDefaultType(64), arg);  \
   M(KmerShapeDefaultType(65), arg);  \
   M(KmerShapeDefaultType(68), arg);  \
   M(KmerShapeDefaultType(80), arg);  \
   M(KmerShapeDefaultType(88), arg);  \
   M(KmerShapeDefaultType(96), arg);  \
   M(KmerShapeDefaultType(100), arg); \
   M(KmerShapeDefaultType(128), arg); \
   M(KmerShapeDefaultType(144), arg); \
   M(KmerShapeDefaultType(192), arg); \
   M(KmerShapeDefaultType(200), arg)\

#define FOR_ALL_KSHAPES_D(M, arg ) \
   M(KmerShapeDefaultType(320), arg); \
   M(KmerShapeDefaultType(368), arg); \
   M(KmerShapeDefaultType(400), arg); \
   M(KmerShapeDefaultType(500), arg); \
   M(KmerShapeDefaultType(544), arg); \
   M(KmerShapeDefaultType(640), arg); \
   M(KmerShapeDefaultType(720), arg); \
   M(KmerShapeDefaultType(1000), arg); \
   M(KmerShapeDefaultType(1200), arg); \
   M(KmerShapeDefaultType(1600), arg); \
   M(KmerShapeDefaultType(2000), arg); \
   M(KmerShapeDefaultType(4000), arg); \
   M(KmerShapeDefaultType(10000), arg)

/* Curently have no use for gapped kmers:

#define FOR_ALL_KSHAPES_E(M, arg) \
   M(KmerShapeMidGapType(8,  1), arg); \
   M(KmerShapeMidGapType(12, 1), arg); \
   M(KmerShapeMidGapType(16, 1), arg); \
   M(KmerShapeMidGapType(20, 1), arg)

#define FOR_ALL_KSHAPES_F(M, arg) \
   M(KmerShapeMidGapType(24, 1), arg); \
   M(KmerShapeMidGapType(26, 1), arg); \
   M(KmerShapeMidGapType(28, 1), arg); \
   M(KmerShapeMidGapType(29, 1), arg); \
   M(KmerShapeMidGapType(32, 1), arg); \
   M(KmerShapeMidGapType(36, 1), arg)

#define FOR_ALL_KSHAPES_G(M, arg) \
   M(KmerShapeMidGapType(8,  4), arg); \
   M(KmerShapeMidGapType(12, 4), arg)

#define FOR_ALL_KSHAPES_H(M, arg) \
   M(KmerShapeMidGapType(16, 4), arg); \
   M(KmerShapeMidGapType(20, 4), arg)

#define FOR_ALL_KSHAPES_I(M, arg) \
   M(KmerShapeMidGapType(24, 4), arg); \
   M(KmerShapeMidGapType(28, 4), arg); \
   M(KmerShapeMidGapType(29, 4), arg)
*/

/**
   Macro: FOR_ALL_K
   
   Call the specified macro for each kmer length that we support.
*/
#define FOR_ALL_K(M, arg) \
   M(4, arg); \
   M(8, arg); \
   M(12, arg); \
   M(16, arg); \
   M(19, arg); \
   M(20, arg); \
   M(21, arg); \
   M(24, arg); \
   M(26, arg); \
   M(28, arg); \
   M(29, arg); \
   M(32, arg); \
   M(36, arg); \
   M(40, arg); \
   M(48, arg); \
   M(49, arg); \
   M(52, arg); \
   M(64, arg); \
   M(65, arg); \
   M(68, arg); \
   M(80, arg); \
   M(88, arg); \
   M(96, arg); \
   M(100, arg); \
   M(128, arg); \
   M(144, arg); \
   M(192, arg); \
   M(200, arg); \
   M(320, arg); \
   M(368, arg); \
   M(400, arg); \
   M(500, arg); \
   M(544, arg); \
   M(640, arg); \
   M(720, arg); \
   M(1000, arg); \
   M(1200, arg); \
   M(1600, arg); \
   M(2000, arg); \
   M(4000, arg); \
   M(10000, arg)

/**
   Macro: FOR_ALL_K_PLUS_1
   
   Call the specified macro for each kmer+1 length that we support.
*/
#define FOR_ALL_K_PLUS_1(M, arg) \
   M(20, arg); \
   M(21, arg); \
   M(29, arg); \
   M(49, arg); \
   M(65, arg)

/**
   Macro: FOR_ALL_K_WITH_K_PLUS_1
   
   Call the specified macro for each kmer length that we support a
   kmer+1 length too.
*/
#define FOR_ALL_K_WITH_K_PLUS_1(M, arg) \
   M(19, arg); \
   M(20, arg); \
   M(28, arg); \
   M(48, arg); \
   M(64, arg)


#else

///////////////////////////
// Put your custom kmer shape set here.  Define just the kmers with which you work.

#define FOR_ALL_KSHAPES_A(M, arg) \
   M(KmerShapeDefaultType(8), arg); \
   M(KmerShapeDefaultType(20), arg); \
   M(KmerShapeDefaultType(24), arg)

#define FOR_ALL_KSHAPES_B(M, arg) \
     M(KmerShapeDefaultType(12), arg); \
     M(KmerShapeDefaultType(16), arg)

#define FOR_ALL_KSHAPES_C(M, arg)    \
   M(KmerShapeMidGapType(16, 1), arg); \
   M(KmerShapeMidGapType(20, 1), arg)

#define FOR_ALL_KSHAPES_D(M, arg) \
   NO_KSHAPES(D) 

/**
   Macro: FOR_ALL_K
   
   Call the specified macro for each kmer length that we support.
*/
#define FOR_ALL_K(M, arg) \
   M(8, arg);  \
   M(12, arg); \
   M(16, arg); \
   M(20, arg); \
   M(21, arg); \
   M(24, arg) 

/**
   Macro: FOR_ALL_K_PLUS_1
   
   Call the specified macro for each kmer+1 length that we support.
*/
#define FOR_ALL_K_PLUS_1(M, arg) \
   M(21, arg)

/**
   Macro: FOR_ALL_K_WITH_K_PLUS_1
   
   Call the specified macro for each kmer length that we support a
   kmer+1 length too.
*/
#define FOR_ALL_K_WITH_K_PLUS_1(M, arg) \
   M(20, arg)


#endif


#endif
// #ifndef __INCLUDE_SupportedKmerShapes_h
