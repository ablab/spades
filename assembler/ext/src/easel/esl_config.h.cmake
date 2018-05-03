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

/* Large file support
 * Must precede any header file inclusion.
 */
#cmakedefine _FILE_OFFSET_BITS
#cmakedefine _LARGE_FILES
#cmakedefine _LARGEFILE_SOURCE

/* Debugging verbosity (0=none;3=most verbose)
 */
#cmakedefine eslDEBUGLEVEL

/* System headers
 */
#cmakedefine HAVE_ENDIAN_H
#cmakedefine HAVE_INTTYPES_H
#cmakedefine HAVE_STDINT_H
#cmakedefine HAVE_UNISTD_H
#cmakedefine HAVE_SYS_TYPES_H
#cmakedefine HAVE_STRINGS_H

#cmakedefine HAVE_SYS_PARAM_H
#cmakedefine HAVE_SYS_SYSCTL_H

#cmakedefine HAVE_EMMINTRIN_H
#cmakedefine HAVE_PMMINTRIN_H
#cmakedefine HAVE_XMMINTRIN_H

#cmakedefine HAVE_ALTIVEC_H

/* Types
 */
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

/* Optional packages
 */
#cmakedefine HAVE_LIBGSL

/* Optional parallel implementation support
 */
#cmakedefine HAVE_SSE2
#cmakedefine HAVE_VMX
#cmakedefine HAVE_MPI
#cmakedefine HAVE_PTHREAD

#cmakedefine HAVE_SSE2_CAST

/* Programs */
#cmakedefine HAVE_GZIP

/* Functions */
#cmakedefine HAVE_CHMOD
#cmakedefine HAVE_FSEEKO
#cmakedefine HAVE_FSTAT
#cmakedefine HAVE_GETCWD
#cmakedefine HAVE_GETPID
#cmakedefine HAVE_MKSTEMP
#cmakedefine HAVE_POPEN
#cmakedefine HAVE_PUTENV
#cmakedefine HAVE_STAT
#cmakedefine HAVE_STRCASECMP
#cmakedefine HAVE_SYSCONF
#cmakedefine HAVE_SYSCTL
#cmakedefine HAVE_TIMES
#cmakedefine HAVE_ERFC

#cmakedefine HAVE_FUNC_ATTRIBUTE_NORETURN // Compiler supports __attribute__ tag, which we use to help w/ clang static analysis.

/* Function behavior */
#define eslSTOPWATCH_HIGHRES

/*****************************************************************
 * Available augmentations.
 *
 * If you grab a single module from Easel to use it by itself,
 * leave all these #undef'd; you have no augmentations.
 *
 * If you grab additional Easel .c files, you can enable any
 * augmentations they provide to other modules by #defining the
 * modules you have below. Alternatively, you can -D them on
 * the compile line, as in cc -DeslAUGMENT_SSI -DeslAUGMENT_MSA.
 *
 * If you compile and install the complete Easel library, all of these
 * get #defined automatically by ./configure, plus the eslLIBRARY flag
 * which means the full library with all augmentations is
 * available. So, if you steal files from an installed library, just
 * set these all back to #undef (depending on which files you have).
 *****************************************************************/
#cmakedefine eslLIBRARY

#ifndef eslLIBRARY
#undef eslAUGMENT_ALPHABET
#undef eslAUGMENT_NCBI
#undef eslAUGMENT_DMATRIX
#undef eslAUGMENT_FILEPARSER
#undef eslAUGMENT_GEV
#undef eslAUGMENT_GUMBEL
#undef eslAUGMENT_HISTOGRAM
#undef eslAUGMENT_KEYHASH
#undef eslAUGMENT_MINIMIZER
#undef eslAUGMENT_MSA
#undef eslAUGMENT_RANDOM
#undef eslAUGMENT_RANDOMSEQ
#undef eslAUGMENT_SSI
#undef eslAUGMENT_STATS
#endif

#ifdef eslLIBRARY
#define eslAUGMENT_ALPHABET
#define eslAUGMENT_NCBI
#define eslAUGMENT_DMATRIX
#define eslAUGMENT_FILEPARSER
#define eslAUGMENT_GEV
#define eslAUGMENT_GUMBEL
#define eslAUGMENT_HISTOGRAM
#define eslAUGMENT_KEYHASH
#define eslAUGMENT_MINIMIZER
#define eslAUGMENT_MSA
#define eslAUGMENT_RANDOM
#define eslAUGMENT_RANDOMSEQ
#define eslAUGMENT_SSI
#define eslAUGMENT_STATS
#endif


#endif /*eslCONFIG_INCLUDED*/
