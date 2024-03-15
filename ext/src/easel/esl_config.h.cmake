/* esl_config.h.in  [input to configure]
 * 
 * System-dependent configuration of Easel, by autoconf.
 * 
 * This file should be included in all Easel .c files before
 * anything else, because it may set #define's that control
 * behaviour of system includes and system libraries. An example
 * is large file support.
 * 
 */
#ifndef eslCONFIG_INCLUDED
#define eslCONFIG_INCLUDED

/* Version info.
 */
#cmakedefine EASEL_VERSION   "@EASEL_VERSION@"
#cmakedefine EASEL_DATE      "@EASEL_DATE@"
#cmakedefine EASEL_COPYRIGHT "@EASEL_COPYRIGHT@"
#cmakedefine EASEL_LICENSE   "@EASEL_LICENSE@"

/* Debugging/assertion hooks & verbosity level (0=none;3=most verbose) */
#cmakedefine eslDEBUGLEVEL

/* Optional parallel implementation support */
#cmakedefine eslENABLE_SSE
#cmakedefine eslENABLE_SSE4
#cmakedefine eslENABLE_AVX
#cmakedefine eslENABLE_AVX512
#cmakedefine eslENABLE_NEON
#cmakedefine eslENABLE_VMX

#cmakedefine eslHAVE_NEON_AARCH64

#cmakedefine eslENABLE_CUDA  // Should we build CUDA acceleration?

#cmakedefine HAVE_FLUSH_ZERO_MODE        // on x86 platforms: we can turn off denormalized floating point math,
#cmakedefine HAVE_DENORMALS_ZERO_MODE    //   which often incurs performance penalty. See simdvec.md in HMMER.

#cmakedefine HAVE_MPI
#cmakedefine HAVE_PTHREAD

/* Programs */
#cmakedefine HAVE_GZIP

/* Libraries */
#cmakedefine HAVE_LIBGSL

/* Headers */
#cmakedefine HAVE_ENDIAN_H
#cmakedefine HAVE_INTTYPES_H
#cmakedefine HAVE_STDINT_H
#cmakedefine HAVE_UNISTD_H
#cmakedefine HAVE_SYS_TYPES_H
#cmakedefine HAVE_STRINGS_H
#cmakedefine HAVE_NETINET_IN_H	/* On FreeBSD, you need netinet/in.h for struct sockaddr_in */

#cmakedefine HAVE_SYS_PARAM_H
#cmakedefine HAVE_SYS_SYSCTL_H

/* Types */
#cmakedefine WORDS_BIGENDIAN
#cmakedefine int8_t
#cmakedefine int16_t
#cmakedefine int32_t
#cmakedefine int64_t
#cmakedefine uint8_t
#cmakedefine uint16_t
#cmakedefine uint32_t
#cmakedefine uint64_t
#cmakedefine off_t

/* Compiler characteristics */
#cmakedefine HAVE_FUNC_ATTRIBUTE_NORETURN // Compiler supports __attribute__((__noreturn__)), helps w/ clang static analysis.
#cmakedefine HAVE_FUNC_ATTRIBUTE_FORMAT   // Compiler supports __attribute__((format(a,b,c))), typechecking printf-like functions

/* Functions */
#cmakedefine HAVE_ALIGNED_ALLOC   // esl_alloc
#cmakedefine HAVE_ERFC            // esl_stats
#cmakedefine HAVE_GETCWD          // esl_getcwd
#cmakedefine HAVE_GETPID          // esl_random
#cmakedefine HAVE__MM_MALLOC      // esl_alloc
#cmakedefine HAVE_POPEN           // various file parsers that check for piped input
#cmakedefine HAVE_POSIX_MEMALIGN  // esl_alloc
#cmakedefine HAVE_STRCASECMP      // easel::esl_strcasecmp()
#cmakedefine HAVE_STRSEP          // easel::esl_strsep()
#cmakedefine HAVE_SYSCONF         // esl_threads, asking system for cpu number
#cmakedefine HAVE_SYSCTL          // esl_threads, ""
#cmakedefine HAVE_TIMES           // esl_stopwatch

#cmakedefine HAVE_FSEEKO

/* System services */
#cmakedefine _FILE_OFFSET_BITS    // Large file support; possibly archaic now?
#cmakedefine _LARGE_FILES         //  ""
#cmakedefine _LARGEFILE_SOURCE    //  ""

 
/* Function behavior */
#define eslSTOPWATCH_HIGHRES

#endif /*eslCONFIG_INCLUDED*/

