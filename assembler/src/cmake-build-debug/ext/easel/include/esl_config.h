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
#define EASEL_VERSION   "0.48"
#define EASEL_DATE      "Nov 2020"
#define EASEL_COPYRIGHT "Copyright (C) 2020 Howard Hughes Medical Institute"
#define EASEL_LICENSE   "Freely distributed under the BSD open source license."

/* Debugging/assertion hooks & verbosity level (0=none;3=most verbose) */
/* #undef eslDEBUGLEVEL */

/* Optional parallel implementation support */
#define eslENABLE_SSE
/* #undef eslENABLE_SSE4 */
/* #undef eslENABLE_AVX */
/* #undef eslENABLE_AVX512 */
/* #undef eslENABLE_NEON */
/* #undef eslENABLE_VMX */

/* #undef eslHAVE_NEON_AARCH64 */

/* #undef eslENABLE_CUDA */

#define HAVE_FLUSH_ZERO_MODE        // on x86 platforms: we can turn off denormalized floating point math,
/* #undef HAVE_DENORMALS_ZERO_MODE */

/* #undef HAVE_MPI */
#define HAVE_PTHREAD

/* Programs */
/* #undef HAVE_GZIP */

/* Libraries */
/* #undef HAVE_LIBGSL */

/* Headers */
#define HAVE_ENDIAN_H
#define HAVE_INTTYPES_H
#define HAVE_STDINT_H
#define HAVE_UNISTD_H
#define HAVE_SYS_TYPES_H
#define HAVE_STRINGS_H
#define HAVE_NETINET_IN_H	/* On FreeBSD, you need netinet/in.h for struct sockaddr_in */

#define HAVE_SYS_PARAM_H
/* #undef HAVE_SYS_SYSCTL_H */

/* Types */
/* #undef WORDS_BIGENDIAN */
/* #undef int8_t */
/* #undef int16_t */
/* #undef int32_t */
/* #undef int64_t */
/* #undef uint8_t */
/* #undef uint16_t */
/* #undef uint32_t */
/* #undef uint64_t */
/* #undef off_t */

/* Compiler characteristics */
#define HAVE_FUNC_ATTRIBUTE_NORETURN // Compiler supports __attribute__((__noreturn__)), helps w/ clang static analysis.
#define HAVE_FUNC_ATTRIBUTE_FORMAT   // Compiler supports __attribute__((format(a,b,c))), typechecking printf-like functions

/* Functions */
#define HAVE_ALIGNED_ALLOC   // esl_alloc
/* #undef HAVE_ERFC */
/* #undef HAVE_GETCWD */
#define HAVE_GETPID          // esl_random
/* #undef HAVE__MM_MALLOC */
#define HAVE_POPEN           // various file parsers that check for piped input
#define HAVE_POSIX_MEMALIGN  // esl_alloc
#define HAVE_STRCASECMP      // easel::esl_strcasecmp()
#define HAVE_STRSEP          // easel::esl_strsep()
#define HAVE_SYSCONF         // esl_threads, asking system for cpu number
/* #undef HAVE_SYSCTL */
#define HAVE_TIMES           // esl_stopwatch

#define HAVE_FSEEKO

/* System services */
/* #undef _FILE_OFFSET_BITS */
/* #undef _LARGE_FILES */
/* #undef _LARGEFILE_SOURCE */

 
/* Function behavior */
#define eslSTOPWATCH_HIGHRES

#endif /*eslCONFIG_INCLUDED*/

