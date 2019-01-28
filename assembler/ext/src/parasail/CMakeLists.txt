CMAKE_MINIMUM_REQUIRED( VERSION 3.2 )

FILE( STRINGS ../../include/parasail/parasail.h MAJOR REGEX "^#define +PARASAIL_VERSION_MAJOR +[0-9]+$" )
FILE( STRINGS ../../include/parasail/parasail.h MINOR REGEX "^#define +PARASAIL_VERSION_MINOR +[0-9]+$" )
FILE( STRINGS ../../include/parasail/parasail.h PATCH REGEX "^#define +PARASAIL_VERSION_PATCH +[0-9]+$" )
STRING( REPLACE "#define PARASAIL_VERSION_MAJOR " "" MAJOR ${MAJOR} )
STRING( REPLACE "#define PARASAIL_VERSION_MINOR " "" MINOR ${MINOR} )
STRING( REPLACE "#define PARASAIL_VERSION_PATCH " "" PATCH ${PATCH} )

PROJECT( parasail VERSION ${MAJOR}.${MINOR}.${PATCH} )

SET( CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake" ${CMAKE_MODULE_PATH} )

#OPTION( BUILD_SHARED_LIBS "Build Shared Libraries" ON )
SET( LIB_TYPE STATIC )
IF( BUILD_SHARED_LIBS )
  SET( LIB_TYPE SHARED )
  # All OBJECT libraries below are later consumed by a shared library so
  # the following is OK. Saves on setting it as a per target property.
  SET( CMAKE_POSITION_INDEPENDENT_CODE ON )
ENDIF( )

INCLUDE( CheckFunctionExists )
INCLUDE( CheckIncludeFile )
INCLUDE( CheckLibraryExists )
INCLUDE( CheckSymbolExists )
INCLUDE( CheckTypeSize )
INCLUDE( CheckCSourceCompiles )
INCLUDE( FindOpenMP )
INCLUDE( FindSSE2 )
INCLUDE( FindSSE41 )
INCLUDE( FindAVX2 )
INCLUDE( Findxgetbv )
INCLUDE( FindAltiVec )
INCLUDE( FindNEON )
INCLUDE( FindZLIB )
INCLUDE( TestBigEndian )

MESSAGE( STATUS "Check if the CPU is POWER" )
IF( ${CMAKE_SYSTEM_PROCESSOR} MATCHES "p(ower)?pc" )
    SET( IS_POWER_ISA TRUE )
    MESSAGE( STATUS "Check if the CPU is POWER - TRUE" )
ELSE( )
    SET( IS_POWER_ISA FALSE )
    MESSAGE( STATUS "Check if the CPU is POWER - FALSE" )
ENDIF( )

MESSAGE( STATUS "Check if the CPU is ARM" )
IF( ${CMAKE_SYSTEM_PROCESSOR} MATCHES "arm.*" OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES "aarch64.*" )
    SET( IS_ARM_ISA TRUE )
    MESSAGE( STATUS "Check if the CPU is ARM - TRUE" )
ELSE( )
    SET( IS_ARM_ISA FALSE )
    MESSAGE( STATUS "Check if the CPU is ARM - FALSE" )
ENDIF( )

CHECK_TYPE_SIZE( "int" SIZEOF_INT )
CHECK_TYPE_SIZE( "int8_t" SIZEOF_INT8_T )
CHECK_TYPE_SIZE( "int16_t" SIZEOF_INT16_T )
CHECK_TYPE_SIZE( "int32_t" SIZEOF_INT32_T )
CHECK_TYPE_SIZE( "int64_t" SIZEOF_INT64_T )
CHECK_TYPE_SIZE( "uint8_t" SIZEOF_UINT8_T )
CHECK_TYPE_SIZE( "uint16_t" SIZEOF_UINT16_T )
CHECK_TYPE_SIZE( "uint32_t" SIZEOF_UINT32_T )
CHECK_TYPE_SIZE( "uint64_t" SIZEOF_UINT64_T )
TEST_BIG_ENDIAN( WORDS_BIGENDIAN )

# Check for restrict keyword
FOREACH( ac_kw __restrict __restrict__ _Restrict restrict )
    CHECK_C_SOURCE_COMPILES(
"
typedef int * int_ptr;
int foo (int_ptr ${ac_kw} ip) {
    return ip[0];
}
int main() {
    int s[1];
    int * ${ac_kw} t = s;
    t[0] = 0;
    return foo(t);
}
"
    HAVE_RESTRICT )
    IF( HAVE_RESTRICT )
        SET( ac_cv_c_restrict ${ac_kw} )
        BREAK( )
    ENDIF( )
ENDFOREACH( )
IF( HAVE_RESTRICT )
    SET( restrict ${ac_cv_c_restrict} )
ELSE( )
    SET( restrict " " )
ENDIF( )

# Check for inline keyword
CHECK_C_SOURCE_COMPILES(
"
typedef int foo_t;
static inline foo_t static_foo(){return 0;}
foo_t foo(){return 0;}
int main(int argc, char *argv[]){return 0;}
"
    HAVE_INLINE_NATIVE )
IF( HAVE_INLINE_NATIVE )
ELSE ( )
    FOREACH( ac_kw __inline__ __inline )
        CHECK_C_SOURCE_COMPILES(
"
typedef int foo_t;
static ${ac_kw} foo_t static_foo(){return 0;}
foo_t foo(){return 0;}
int main(int argc, char *argv[]){return 0;}
"
        HAVE_INLINE )
        IF( HAVE_INLINE )
            SET( ac_cv_c_inline ${ac_kw} )
            BREAK( )
        ENDIF( )
    ENDFOREACH( )
    IF( HAVE_INLINE )
        SET( inline ${ac_cv_c_inline} )
    ELSE( )
        SET( inline " " )
    ENDIF( )
ENDIF( )

SET( CMAKE_REQUIRED_INCLUDES stdint.h )

INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_LIST_DIR} )
INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_LIST_DIR}/contrib )
INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_LIST_DIR}/contrib/sais-lite-lcp )
INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_LIST_DIR}/contrib/ssw )
INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_LIST_DIR}/src )
INCLUDE_DIRECTORIES( ${CMAKE_BINARY_DIR} )
IF( NEON_FOUND )
INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_LIST_DIR}/contrib/simde/simde )
ENDIF( )

CHECK_INCLUDE_FILE( sys/mman.h HAVE_SYS_MMAN_H )
CHECK_INCLUDE_FILE( unistd.h HAVE_UNISTD_H )
CHECK_INCLUDE_FILE( windows.h HAVE_WINDOWS_H )

#CHECK_FUNCTION_EXISTS( _aligned_malloc HAVE__ALIGNED_MALLOC )
#CHECK_FUNCTION_EXISTS( posix_memalign HAVE_POSIX_MEMALIGN )
#CHECK_FUNCTION_EXISTS( aligned_alloc HAVE_ALIGNED_ALLOC )
#CHECK_FUNCTION_EXISTS( memalign HAVE_MEMALIGN )
CHECK_SYMBOL_EXISTS( _aligned_malloc malloc.h HAVE__ALIGNED_MALLOC )
CHECK_SYMBOL_EXISTS( posix_memalign stdlib.h HAVE_POSIX_MEMALIGN )
CHECK_SYMBOL_EXISTS( aligned_alloc stdlib.h HAVE_ALIGNED_ALLOC )
CHECK_SYMBOL_EXISTS( memalign malloc.h HAVE_MEMALIGN )
CHECK_FUNCTION_EXISTS( getopt HAVE_GETOPT )
CHECK_FUNCTION_EXISTS( poll HAVE_POLL )
#CHECK_FUNCTION_EXISTS( filelength HAVE_FILELENGTH )
#CHECK_FUNCTION_EXISTS( GetSystemTimeAsFileTime HAVE_GETSYSTEMTIMEASFILETIME )
#CHECK_FUNCTION_EXISTS( SetUnhandledExceptionFilter HAVE_SETUNHANDLEDEXCEPTIONFILTER )
CHECK_SYMBOL_EXISTS( filelength io.h HAVE_FILELENGTH )
CHECK_SYMBOL_EXISTS( GetSystemTimeAsFileTime windows.h HAVE_GETSYSTEMTIMEASFILETIME )
CHECK_SYMBOL_EXISTS( SetUnhandledExceptionFilter windows.h HAVE_SETUNHANDLEDEXCEPTIONFILTER )
CHECK_FUNCTION_EXISTS( clock_get_time HAVE_CLOCK_GET_TIME_NATIVE )
CHECK_FUNCTION_EXISTS( clock_gettime HAVE_CLOCK_GETTIME_NATIVE )
CHECK_LIBRARY_EXISTS( rt clock_gettime "" HAVE_CLOCK_GETTIME_RT )
CHECK_LIBRARY_EXISTS( pthread clock_gettime "" HAVE_CLOCK_GETTIME_PTHREAD )
CHECK_FUNCTION_EXISTS( sqrt HAVE_SQRT_NATIVE )
CHECK_LIBRARY_EXISTS( m sqrt "" HAVE_SQRT_M )

CONFIGURE_FILE( ${CMAKE_CURRENT_LIST_DIR}/cmake/config.h.in ${CMAKE_BINARY_DIR}/config.h )

SET( SRC_CORE
    src/cigar.c
    src/function_lookup.c
    src/io.c
    src/isastubs.c
    src/matrix_lookup.c
    src/memory.c
    src/parser.c
    src/pssw.c
    src/time.c
    src/nw_dispatch.c
    src/sg_dispatch.c
    src/sw_dispatch.c
    src/dispatch_profile.c
    src/satcheck.c
    src/striped_unwind.c
    src/traceback.c
)

IF( IS_POWER_ISA )
    SET( SRC_CORE ${SRC_CORE} src/cpuid_ppc.c )
ELSEIF( IS_ARM_ISA )
    SET( SRC_CORE ${SRC_CORE} src/cpuid_arm.c )
ELSE( )
    SET( SRC_CORE ${SRC_CORE} src/cpuid.c )
ENDIF( )

SET( SRC_NOVEC
    src/nw.c
    src/sg.c
    src/sw.c
    src/nw_banded.c
    src/nw_scan.c
    src/sg_scan.c
    src/sw_scan.c
    src/nw_stats.c
    src/sg_stats.c
    src/sw_stats.c
    src/nw_stats_scan.c
    src/sg_stats_scan.c
    src/sw_stats_scan.c
)

SET( SRC_TRACE_NOVEC
    src/nw_trace.c
    src/sg_trace.c
    src/sw_trace.c
    src/nw_trace_scan.c
    src/sg_trace_scan.c
    src/sw_trace_scan.c
)

SET( SRC_SSE2
    src/nw_scan_sse2_128_64.c
    src/sg_scan_sse2_128_64.c
    src/sw_scan_sse2_128_64.c
    src/nw_scan_sse2_128_32.c
    src/sg_scan_sse2_128_32.c
    src/sw_scan_sse2_128_32.c
    src/nw_scan_sse2_128_16.c
    src/sg_scan_sse2_128_16.c
    src/sw_scan_sse2_128_16.c
    src/nw_scan_sse2_128_8.c
    src/sg_scan_sse2_128_8.c
    src/sw_scan_sse2_128_8.c
    src/nw_diag_sse2_128_64.c
    src/sg_diag_sse2_128_64.c
    src/sw_diag_sse2_128_64.c
    src/nw_diag_sse2_128_32.c
    src/sg_diag_sse2_128_32.c
    src/sw_diag_sse2_128_32.c
    src/nw_diag_sse2_128_16.c
    src/sg_diag_sse2_128_16.c
    src/sw_diag_sse2_128_16.c
    src/nw_diag_sse2_128_8.c
    src/sg_diag_sse2_128_8.c
    src/sw_diag_sse2_128_8.c
    src/nw_striped_sse2_128_64.c
    src/sg_striped_sse2_128_64.c
    src/sw_striped_sse2_128_64.c
    src/nw_striped_sse2_128_32.c
    src/sg_striped_sse2_128_32.c
    src/sw_striped_sse2_128_32.c
    src/nw_striped_sse2_128_16.c
    src/sg_striped_sse2_128_16.c
    src/sw_striped_sse2_128_16.c
    src/nw_striped_sse2_128_8.c
    src/sg_striped_sse2_128_8.c
    src/sw_striped_sse2_128_8.c
    src/nw_stats_scan_sse2_128_64.c
    src/sg_stats_scan_sse2_128_64.c
    src/sw_stats_scan_sse2_128_64.c
    src/nw_stats_scan_sse2_128_32.c
    src/sg_stats_scan_sse2_128_32.c
    src/sw_stats_scan_sse2_128_32.c
    src/nw_stats_scan_sse2_128_16.c
    src/sg_stats_scan_sse2_128_16.c
    src/sw_stats_scan_sse2_128_16.c
    src/nw_stats_scan_sse2_128_8.c
    src/sg_stats_scan_sse2_128_8.c
    src/sw_stats_scan_sse2_128_8.c
    src/nw_stats_diag_sse2_128_64.c
    src/sg_stats_diag_sse2_128_64.c
    src/sw_stats_diag_sse2_128_64.c
    src/nw_stats_diag_sse2_128_32.c
    src/sg_stats_diag_sse2_128_32.c
    src/sw_stats_diag_sse2_128_32.c
    src/nw_stats_diag_sse2_128_16.c
    src/sg_stats_diag_sse2_128_16.c
    src/sw_stats_diag_sse2_128_16.c
    src/nw_stats_diag_sse2_128_8.c
    src/sg_stats_diag_sse2_128_8.c
    src/sw_stats_diag_sse2_128_8.c
    src/nw_stats_striped_sse2_128_64.c
    src/sg_stats_striped_sse2_128_64.c
    src/sw_stats_striped_sse2_128_64.c
    src/nw_stats_striped_sse2_128_32.c
    src/sg_stats_striped_sse2_128_32.c
    src/sw_stats_striped_sse2_128_32.c
    src/nw_stats_striped_sse2_128_16.c
    src/sg_stats_striped_sse2_128_16.c
    src/sw_stats_striped_sse2_128_16.c
    src/nw_stats_striped_sse2_128_8.c
    src/sg_stats_striped_sse2_128_8.c
    src/sw_stats_striped_sse2_128_8.c
)

SET( SRC_TRACE_SSE2
    src/nw_trace_scan_sse2_128_64.c
    src/sg_trace_scan_sse2_128_64.c
    src/sw_trace_scan_sse2_128_64.c
    src/nw_trace_scan_sse2_128_32.c
    src/sg_trace_scan_sse2_128_32.c
    src/sw_trace_scan_sse2_128_32.c
    src/nw_trace_scan_sse2_128_16.c
    src/sg_trace_scan_sse2_128_16.c
    src/sw_trace_scan_sse2_128_16.c
    src/nw_trace_scan_sse2_128_8.c
    src/sg_trace_scan_sse2_128_8.c
    src/sw_trace_scan_sse2_128_8.c
    src/nw_trace_diag_sse2_128_64.c
    src/sg_trace_diag_sse2_128_64.c
    src/sw_trace_diag_sse2_128_64.c
    src/nw_trace_diag_sse2_128_32.c
    src/sg_trace_diag_sse2_128_32.c
    src/sw_trace_diag_sse2_128_32.c
    src/nw_trace_diag_sse2_128_16.c
    src/sg_trace_diag_sse2_128_16.c
    src/sw_trace_diag_sse2_128_16.c
    src/nw_trace_diag_sse2_128_8.c
    src/sg_trace_diag_sse2_128_8.c
    src/sw_trace_diag_sse2_128_8.c
    src/nw_trace_striped_sse2_128_64.c
    src/sg_trace_striped_sse2_128_64.c
    src/sw_trace_striped_sse2_128_64.c
    src/nw_trace_striped_sse2_128_32.c
    src/sg_trace_striped_sse2_128_32.c
    src/sw_trace_striped_sse2_128_32.c
    src/nw_trace_striped_sse2_128_16.c
    src/sg_trace_striped_sse2_128_16.c
    src/sw_trace_striped_sse2_128_16.c
    src/nw_trace_striped_sse2_128_8.c
    src/sg_trace_striped_sse2_128_8.c
    src/sw_trace_striped_sse2_128_8.c
)

SET( SRC_SSE41
    src/nw_scan_sse41_128_64.c
    src/sg_scan_sse41_128_64.c
    src/sw_scan_sse41_128_64.c
    src/nw_scan_sse41_128_32.c
    src/sg_scan_sse41_128_32.c
    src/sw_scan_sse41_128_32.c
    src/nw_scan_sse41_128_16.c
    src/sg_scan_sse41_128_16.c
    src/sw_scan_sse41_128_16.c
    src/nw_scan_sse41_128_8.c
    src/sg_scan_sse41_128_8.c
    src/sw_scan_sse41_128_8.c
    src/nw_diag_sse41_128_64.c
    src/sg_diag_sse41_128_64.c
    src/sw_diag_sse41_128_64.c
    src/nw_diag_sse41_128_32.c
    src/sg_diag_sse41_128_32.c
    src/sw_diag_sse41_128_32.c
    src/nw_diag_sse41_128_16.c
    src/sg_diag_sse41_128_16.c
    src/sw_diag_sse41_128_16.c
    src/nw_diag_sse41_128_8.c
    src/sg_diag_sse41_128_8.c
    src/sw_diag_sse41_128_8.c
    src/nw_striped_sse41_128_64.c
    src/sg_striped_sse41_128_64.c
    src/sw_striped_sse41_128_64.c
    src/nw_striped_sse41_128_32.c
    src/sg_striped_sse41_128_32.c
    src/sw_striped_sse41_128_32.c
    src/nw_striped_sse41_128_16.c
    src/sg_striped_sse41_128_16.c
    src/sw_striped_sse41_128_16.c
    src/nw_striped_sse41_128_8.c
    src/sg_striped_sse41_128_8.c
    src/sw_striped_sse41_128_8.c
    src/sw_blocked_sse41_128_32.c
    src/sw_blocked_sse41_128_16.c
    src/nw_stats_scan_sse41_128_64.c
    src/sg_stats_scan_sse41_128_64.c
    src/sw_stats_scan_sse41_128_64.c
    src/nw_stats_scan_sse41_128_32.c
    src/sg_stats_scan_sse41_128_32.c
    src/sw_stats_scan_sse41_128_32.c
    src/nw_stats_scan_sse41_128_16.c
    src/sg_stats_scan_sse41_128_16.c
    src/sw_stats_scan_sse41_128_16.c
    src/nw_stats_scan_sse41_128_8.c
    src/sg_stats_scan_sse41_128_8.c
    src/sw_stats_scan_sse41_128_8.c
    src/nw_stats_diag_sse41_128_64.c
    src/sg_stats_diag_sse41_128_64.c
    src/sw_stats_diag_sse41_128_64.c
    src/nw_stats_diag_sse41_128_32.c
    src/sg_stats_diag_sse41_128_32.c
    src/sw_stats_diag_sse41_128_32.c
    src/nw_stats_diag_sse41_128_16.c
    src/sg_stats_diag_sse41_128_16.c
    src/sw_stats_diag_sse41_128_16.c
    src/nw_stats_diag_sse41_128_8.c
    src/sg_stats_diag_sse41_128_8.c
    src/sw_stats_diag_sse41_128_8.c
    src/nw_stats_striped_sse41_128_64.c
    src/sg_stats_striped_sse41_128_64.c
    src/sw_stats_striped_sse41_128_64.c
    src/nw_stats_striped_sse41_128_32.c
    src/sg_stats_striped_sse41_128_32.c
    src/sw_stats_striped_sse41_128_32.c
    src/nw_stats_striped_sse41_128_16.c
    src/sg_stats_striped_sse41_128_16.c
    src/sw_stats_striped_sse41_128_16.c
    src/nw_stats_striped_sse41_128_8.c
    src/sg_stats_striped_sse41_128_8.c
    src/sw_stats_striped_sse41_128_8.c
)

SET( SRC_TRACE_SSE41
    src/nw_trace_scan_sse41_128_64.c
    src/sg_trace_scan_sse41_128_64.c
    src/sw_trace_scan_sse41_128_64.c
    src/nw_trace_scan_sse41_128_32.c
    src/sg_trace_scan_sse41_128_32.c
    src/sw_trace_scan_sse41_128_32.c
    src/nw_trace_scan_sse41_128_16.c
    src/sg_trace_scan_sse41_128_16.c
    src/sw_trace_scan_sse41_128_16.c
    src/nw_trace_scan_sse41_128_8.c
    src/sg_trace_scan_sse41_128_8.c
    src/sw_trace_scan_sse41_128_8.c
    src/nw_trace_diag_sse41_128_64.c
    src/sg_trace_diag_sse41_128_64.c
    src/sw_trace_diag_sse41_128_64.c
    src/nw_trace_diag_sse41_128_32.c
    src/sg_trace_diag_sse41_128_32.c
    src/sw_trace_diag_sse41_128_32.c
    src/nw_trace_diag_sse41_128_16.c
    src/sg_trace_diag_sse41_128_16.c
    src/sw_trace_diag_sse41_128_16.c
    src/nw_trace_diag_sse41_128_8.c
    src/sg_trace_diag_sse41_128_8.c
    src/sw_trace_diag_sse41_128_8.c
    src/nw_trace_striped_sse41_128_64.c
    src/sg_trace_striped_sse41_128_64.c
    src/sw_trace_striped_sse41_128_64.c
    src/nw_trace_striped_sse41_128_32.c
    src/sg_trace_striped_sse41_128_32.c
    src/sw_trace_striped_sse41_128_32.c
    src/nw_trace_striped_sse41_128_16.c
    src/sg_trace_striped_sse41_128_16.c
    src/sw_trace_striped_sse41_128_16.c
    src/nw_trace_striped_sse41_128_8.c
    src/sg_trace_striped_sse41_128_8.c
    src/sw_trace_striped_sse41_128_8.c
)

SET( SRC_AVX2
    src/nw_scan_avx2_256_64.c
    src/sg_scan_avx2_256_64.c
    src/sw_scan_avx2_256_64.c
    src/nw_scan_avx2_256_32.c
    src/sg_scan_avx2_256_32.c
    src/sw_scan_avx2_256_32.c
    src/nw_scan_avx2_256_16.c
    src/sg_scan_avx2_256_16.c
    src/sw_scan_avx2_256_16.c
    src/nw_scan_avx2_256_8.c
    src/sg_scan_avx2_256_8.c
    src/sw_scan_avx2_256_8.c
    src/nw_diag_avx2_256_64.c
    src/sg_diag_avx2_256_64.c
    src/sw_diag_avx2_256_64.c
    src/nw_diag_avx2_256_32.c
    src/sg_diag_avx2_256_32.c
    src/sw_diag_avx2_256_32.c
    src/nw_diag_avx2_256_16.c
    src/sg_diag_avx2_256_16.c
    src/sw_diag_avx2_256_16.c
    src/nw_diag_avx2_256_8.c
    src/sg_diag_avx2_256_8.c
    src/sw_diag_avx2_256_8.c
    src/nw_striped_avx2_256_64.c
    src/sg_striped_avx2_256_64.c
    src/sw_striped_avx2_256_64.c
    src/nw_striped_avx2_256_32.c
    src/sg_striped_avx2_256_32.c
    src/sw_striped_avx2_256_32.c
    src/nw_striped_avx2_256_16.c
    src/sg_striped_avx2_256_16.c
    src/sw_striped_avx2_256_16.c
    src/nw_striped_avx2_256_8.c
    src/sg_striped_avx2_256_8.c
    src/sw_striped_avx2_256_8.c
    src/nw_stats_scan_avx2_256_64.c
    src/sg_stats_scan_avx2_256_64.c
    src/sw_stats_scan_avx2_256_64.c
    src/nw_stats_scan_avx2_256_32.c
    src/sg_stats_scan_avx2_256_32.c
    src/sw_stats_scan_avx2_256_32.c
    src/nw_stats_scan_avx2_256_16.c
    src/sg_stats_scan_avx2_256_16.c
    src/sw_stats_scan_avx2_256_16.c
    src/nw_stats_scan_avx2_256_8.c
    src/sg_stats_scan_avx2_256_8.c
    src/sw_stats_scan_avx2_256_8.c
    src/nw_stats_diag_avx2_256_64.c
    src/sg_stats_diag_avx2_256_64.c
    src/sw_stats_diag_avx2_256_64.c
    src/nw_stats_diag_avx2_256_32.c
    src/sg_stats_diag_avx2_256_32.c
    src/sw_stats_diag_avx2_256_32.c
    src/nw_stats_diag_avx2_256_16.c
    src/sg_stats_diag_avx2_256_16.c
    src/sw_stats_diag_avx2_256_16.c
    src/nw_stats_diag_avx2_256_8.c
    src/sg_stats_diag_avx2_256_8.c
    src/sw_stats_diag_avx2_256_8.c
    src/nw_stats_striped_avx2_256_64.c
    src/sg_stats_striped_avx2_256_64.c
    src/sw_stats_striped_avx2_256_64.c
    src/nw_stats_striped_avx2_256_32.c
    src/sg_stats_striped_avx2_256_32.c
    src/sw_stats_striped_avx2_256_32.c
    src/nw_stats_striped_avx2_256_16.c
    src/sg_stats_striped_avx2_256_16.c
    src/sw_stats_striped_avx2_256_16.c
    src/nw_stats_striped_avx2_256_8.c
    src/sg_stats_striped_avx2_256_8.c
    src/sw_stats_striped_avx2_256_8.c
)

SET( SRC_TRACE_AVX2
    src/nw_trace_scan_avx2_256_64.c
    src/sg_trace_scan_avx2_256_64.c
    src/sw_trace_scan_avx2_256_64.c
    src/nw_trace_scan_avx2_256_32.c
    src/sg_trace_scan_avx2_256_32.c
    src/sw_trace_scan_avx2_256_32.c
    src/nw_trace_scan_avx2_256_16.c
    src/sg_trace_scan_avx2_256_16.c
    src/sw_trace_scan_avx2_256_16.c
    src/nw_trace_scan_avx2_256_8.c
    src/sg_trace_scan_avx2_256_8.c
    src/sw_trace_scan_avx2_256_8.c
    src/nw_trace_diag_avx2_256_64.c
    src/sg_trace_diag_avx2_256_64.c
    src/sw_trace_diag_avx2_256_64.c
    src/nw_trace_diag_avx2_256_32.c
    src/sg_trace_diag_avx2_256_32.c
    src/sw_trace_diag_avx2_256_32.c
    src/nw_trace_diag_avx2_256_16.c
    src/sg_trace_diag_avx2_256_16.c
    src/sw_trace_diag_avx2_256_16.c
    src/nw_trace_diag_avx2_256_8.c
    src/sg_trace_diag_avx2_256_8.c
    src/sw_trace_diag_avx2_256_8.c
    src/nw_trace_striped_avx2_256_64.c
    src/sg_trace_striped_avx2_256_64.c
    src/sw_trace_striped_avx2_256_64.c
    src/nw_trace_striped_avx2_256_32.c
    src/sg_trace_striped_avx2_256_32.c
    src/sw_trace_striped_avx2_256_32.c
    src/nw_trace_striped_avx2_256_16.c
    src/sg_trace_striped_avx2_256_16.c
    src/sw_trace_striped_avx2_256_16.c
    src/nw_trace_striped_avx2_256_8.c
    src/sg_trace_striped_avx2_256_8.c
    src/sw_trace_striped_avx2_256_8.c
)

SET( SRC_ALTIVEC
    src/nw_scan_altivec_128_64.c
    src/sg_scan_altivec_128_64.c
    src/sw_scan_altivec_128_64.c
    src/nw_scan_altivec_128_32.c
    src/sg_scan_altivec_128_32.c
    src/sw_scan_altivec_128_32.c
    src/nw_scan_altivec_128_16.c
    src/sg_scan_altivec_128_16.c
    src/sw_scan_altivec_128_16.c
    src/nw_scan_altivec_128_8.c
    src/sg_scan_altivec_128_8.c
    src/sw_scan_altivec_128_8.c
    src/nw_diag_altivec_128_64.c
    src/sg_diag_altivec_128_64.c
    src/sw_diag_altivec_128_64.c
    src/nw_diag_altivec_128_32.c
    src/sg_diag_altivec_128_32.c
    src/sw_diag_altivec_128_32.c
    src/nw_diag_altivec_128_16.c
    src/sg_diag_altivec_128_16.c
    src/sw_diag_altivec_128_16.c
    src/nw_diag_altivec_128_8.c
    src/sg_diag_altivec_128_8.c
    src/sw_diag_altivec_128_8.c
    src/nw_striped_altivec_128_64.c
    src/sg_striped_altivec_128_64.c
    src/sw_striped_altivec_128_64.c
    src/nw_striped_altivec_128_32.c
    src/sg_striped_altivec_128_32.c
    src/sw_striped_altivec_128_32.c
    src/nw_striped_altivec_128_16.c
    src/sg_striped_altivec_128_16.c
    src/sw_striped_altivec_128_16.c
    src/nw_striped_altivec_128_8.c
    src/sg_striped_altivec_128_8.c
    src/sw_striped_altivec_128_8.c
    src/nw_stats_scan_altivec_128_64.c
    src/sg_stats_scan_altivec_128_64.c
    src/sw_stats_scan_altivec_128_64.c
    src/nw_stats_scan_altivec_128_32.c
    src/sg_stats_scan_altivec_128_32.c
    src/sw_stats_scan_altivec_128_32.c
    src/nw_stats_scan_altivec_128_16.c
    src/sg_stats_scan_altivec_128_16.c
    src/sw_stats_scan_altivec_128_16.c
    src/nw_stats_scan_altivec_128_8.c
    src/sg_stats_scan_altivec_128_8.c
    src/sw_stats_scan_altivec_128_8.c
    src/nw_stats_diag_altivec_128_64.c
    src/sg_stats_diag_altivec_128_64.c
    src/sw_stats_diag_altivec_128_64.c
    src/nw_stats_diag_altivec_128_32.c
    src/sg_stats_diag_altivec_128_32.c
    src/sw_stats_diag_altivec_128_32.c
    src/nw_stats_diag_altivec_128_16.c
    src/sg_stats_diag_altivec_128_16.c
    src/sw_stats_diag_altivec_128_16.c
    src/nw_stats_diag_altivec_128_8.c
    src/sg_stats_diag_altivec_128_8.c
    src/sw_stats_diag_altivec_128_8.c
    src/nw_stats_striped_altivec_128_64.c
    src/sg_stats_striped_altivec_128_64.c
    src/sw_stats_striped_altivec_128_64.c
    src/nw_stats_striped_altivec_128_32.c
    src/sg_stats_striped_altivec_128_32.c
    src/sw_stats_striped_altivec_128_32.c
    src/nw_stats_striped_altivec_128_16.c
    src/sg_stats_striped_altivec_128_16.c
    src/sw_stats_striped_altivec_128_16.c
    src/nw_stats_striped_altivec_128_8.c
    src/sg_stats_striped_altivec_128_8.c
    src/sw_stats_striped_altivec_128_8.c
)

SET( SRC_TRACE_ALTIVEC
    src/nw_trace_scan_altivec_128_64.c
    src/sg_trace_scan_altivec_128_64.c
    src/sw_trace_scan_altivec_128_64.c
    src/nw_trace_scan_altivec_128_32.c
    src/sg_trace_scan_altivec_128_32.c
    src/sw_trace_scan_altivec_128_32.c
    src/nw_trace_scan_altivec_128_16.c
    src/sg_trace_scan_altivec_128_16.c
    src/sw_trace_scan_altivec_128_16.c
    src/nw_trace_scan_altivec_128_8.c
    src/sg_trace_scan_altivec_128_8.c
    src/sw_trace_scan_altivec_128_8.c
    src/nw_trace_diag_altivec_128_64.c
    src/sg_trace_diag_altivec_128_64.c
    src/sw_trace_diag_altivec_128_64.c
    src/nw_trace_diag_altivec_128_32.c
    src/sg_trace_diag_altivec_128_32.c
    src/sw_trace_diag_altivec_128_32.c
    src/nw_trace_diag_altivec_128_16.c
    src/sg_trace_diag_altivec_128_16.c
    src/sw_trace_diag_altivec_128_16.c
    src/nw_trace_diag_altivec_128_8.c
    src/sg_trace_diag_altivec_128_8.c
    src/sw_trace_diag_altivec_128_8.c
    src/nw_trace_striped_altivec_128_64.c
    src/sg_trace_striped_altivec_128_64.c
    src/sw_trace_striped_altivec_128_64.c
    src/nw_trace_striped_altivec_128_32.c
    src/sg_trace_striped_altivec_128_32.c
    src/sw_trace_striped_altivec_128_32.c
    src/nw_trace_striped_altivec_128_16.c
    src/sg_trace_striped_altivec_128_16.c
    src/sw_trace_striped_altivec_128_16.c
    src/nw_trace_striped_altivec_128_8.c
    src/sg_trace_striped_altivec_128_8.c
    src/sw_trace_striped_altivec_128_8.c
)

SET( SRC_NEON
    src/nw_scan_neon_128_64.c
    src/sg_scan_neon_128_64.c
    src/sw_scan_neon_128_64.c
    src/nw_scan_neon_128_32.c
    src/sg_scan_neon_128_32.c
    src/sw_scan_neon_128_32.c
    src/nw_scan_neon_128_16.c
    src/sg_scan_neon_128_16.c
    src/sw_scan_neon_128_16.c
    src/nw_scan_neon_128_8.c
    src/sg_scan_neon_128_8.c
    src/sw_scan_neon_128_8.c
    src/nw_diag_neon_128_64.c
    src/sg_diag_neon_128_64.c
    src/sw_diag_neon_128_64.c
    src/nw_diag_neon_128_32.c
    src/sg_diag_neon_128_32.c
    src/sw_diag_neon_128_32.c
    src/nw_diag_neon_128_16.c
    src/sg_diag_neon_128_16.c
    src/sw_diag_neon_128_16.c
    src/nw_diag_neon_128_8.c
    src/sg_diag_neon_128_8.c
    src/sw_diag_neon_128_8.c
    src/nw_striped_neon_128_64.c
    src/sg_striped_neon_128_64.c
    src/sw_striped_neon_128_64.c
    src/nw_striped_neon_128_32.c
    src/sg_striped_neon_128_32.c
    src/sw_striped_neon_128_32.c
    src/nw_striped_neon_128_16.c
    src/sg_striped_neon_128_16.c
    src/sw_striped_neon_128_16.c
    src/nw_striped_neon_128_8.c
    src/sg_striped_neon_128_8.c
    src/sw_striped_neon_128_8.c
    src/nw_stats_scan_neon_128_64.c
    src/sg_stats_scan_neon_128_64.c
    src/sw_stats_scan_neon_128_64.c
    src/nw_stats_scan_neon_128_32.c
    src/sg_stats_scan_neon_128_32.c
    src/sw_stats_scan_neon_128_32.c
    src/nw_stats_scan_neon_128_16.c
    src/sg_stats_scan_neon_128_16.c
    src/sw_stats_scan_neon_128_16.c
    src/nw_stats_scan_neon_128_8.c
    src/sg_stats_scan_neon_128_8.c
    src/sw_stats_scan_neon_128_8.c
    src/nw_stats_diag_neon_128_64.c
    src/sg_stats_diag_neon_128_64.c
    src/sw_stats_diag_neon_128_64.c
    src/nw_stats_diag_neon_128_32.c
    src/sg_stats_diag_neon_128_32.c
    src/sw_stats_diag_neon_128_32.c
    src/nw_stats_diag_neon_128_16.c
    src/sg_stats_diag_neon_128_16.c
    src/sw_stats_diag_neon_128_16.c
    src/nw_stats_diag_neon_128_8.c
    src/sg_stats_diag_neon_128_8.c
    src/sw_stats_diag_neon_128_8.c
    src/nw_stats_striped_neon_128_64.c
    src/sg_stats_striped_neon_128_64.c
    src/sw_stats_striped_neon_128_64.c
    src/nw_stats_striped_neon_128_32.c
    src/sg_stats_striped_neon_128_32.c
    src/sw_stats_striped_neon_128_32.c
    src/nw_stats_striped_neon_128_16.c
    src/sg_stats_striped_neon_128_16.c
    src/sw_stats_striped_neon_128_16.c
    src/nw_stats_striped_neon_128_8.c
    src/sg_stats_striped_neon_128_8.c
    src/sw_stats_striped_neon_128_8.c
)

SET( SRC_TRACE_NEON
    src/nw_trace_scan_neon_128_64.c
    src/sg_trace_scan_neon_128_64.c
    src/sw_trace_scan_neon_128_64.c
    src/nw_trace_scan_neon_128_32.c
    src/sg_trace_scan_neon_128_32.c
    src/sw_trace_scan_neon_128_32.c
    src/nw_trace_scan_neon_128_16.c
    src/sg_trace_scan_neon_128_16.c
    src/sw_trace_scan_neon_128_16.c
    src/nw_trace_scan_neon_128_8.c
    src/sg_trace_scan_neon_128_8.c
    src/sw_trace_scan_neon_128_8.c
    src/nw_trace_diag_neon_128_64.c
    src/sg_trace_diag_neon_128_64.c
    src/sw_trace_diag_neon_128_64.c
    src/nw_trace_diag_neon_128_32.c
    src/sg_trace_diag_neon_128_32.c
    src/sw_trace_diag_neon_128_32.c
    src/nw_trace_diag_neon_128_16.c
    src/sg_trace_diag_neon_128_16.c
    src/sw_trace_diag_neon_128_16.c
    src/nw_trace_diag_neon_128_8.c
    src/sg_trace_diag_neon_128_8.c
    src/sw_trace_diag_neon_128_8.c
    src/nw_trace_striped_neon_128_64.c
    src/sg_trace_striped_neon_128_64.c
    src/sw_trace_striped_neon_128_64.c
    src/nw_trace_striped_neon_128_32.c
    src/sg_trace_striped_neon_128_32.c
    src/sw_trace_striped_neon_128_32.c
    src/nw_trace_striped_neon_128_16.c
    src/sg_trace_striped_neon_128_16.c
    src/sw_trace_striped_neon_128_16.c
    src/nw_trace_striped_neon_128_8.c
    src/sg_trace_striped_neon_128_8.c
    src/sw_trace_striped_neon_128_8.c
)

ADD_LIBRARY( parasail_core OBJECT ${SRC_CORE} )

ADD_LIBRARY( parasail_novec OBJECT ${SRC_NOVEC} )
ADD_LIBRARY( parasail_novec_table OBJECT ${SRC_NOVEC} )
ADD_LIBRARY( parasail_novec_rowcol OBJECT ${SRC_NOVEC} )
ADD_LIBRARY( parasail_novec_trace OBJECT ${SRC_TRACE_NOVEC} )
SET_TARGET_PROPERTIES( parasail_novec_table PROPERTIES COMPILE_DEFINITIONS PARASAIL_TABLE )
SET_TARGET_PROPERTIES( parasail_novec_rowcol PROPERTIES COMPILE_DEFINITIONS PARASAIL_ROWCOL )

IF( SSE2_FOUND )
    ADD_LIBRARY( parasail_sse2 OBJECT ${SRC_SSE2} src/memory_sse.c )
    ADD_LIBRARY( parasail_sse2_table OBJECT ${SRC_SSE2} )
    ADD_LIBRARY( parasail_sse2_rowcol OBJECT ${SRC_SSE2} )
    ADD_LIBRARY( parasail_sse2_trace OBJECT ${SRC_TRACE_SSE2} )
    SET_TARGET_PROPERTIES( parasail_sse2 PROPERTIES COMPILE_FLAGS ${SSE2_C_FLAGS} )
    SET_TARGET_PROPERTIES( parasail_sse2_table PROPERTIES COMPILE_FLAGS ${SSE2_C_FLAGS} )
    SET_TARGET_PROPERTIES( parasail_sse2_rowcol PROPERTIES COMPILE_FLAGS ${SSE2_C_FLAGS} )
    SET_TARGET_PROPERTIES( parasail_sse2_trace PROPERTIES COMPILE_FLAGS ${SSE2_C_FLAGS} )
ELSE( )
    ADD_LIBRARY( parasail_sse2 OBJECT cmake/sse2_dummy.c )
    ADD_LIBRARY( parasail_sse2_table OBJECT cmake/sse2_dummy.c )
    ADD_LIBRARY( parasail_sse2_rowcol OBJECT cmake/sse2_dummy.c )
    ADD_LIBRARY( parasail_sse2_trace OBJECT cmake/sse2_dummy.c )
ENDIF( )
SET_TARGET_PROPERTIES( parasail_sse2_table PROPERTIES COMPILE_DEFINITIONS PARASAIL_TABLE )
SET_TARGET_PROPERTIES( parasail_sse2_rowcol PROPERTIES COMPILE_DEFINITIONS PARASAIL_ROWCOL )
SET_TARGET_PROPERTIES( parasail_sse2_trace PROPERTIES COMPILE_DEFINITIONS PARASAIL_TRACE )

IF( SSE41_FOUND )
    ADD_LIBRARY( parasail_sse41 OBJECT ${SRC_SSE41} )
    ADD_LIBRARY( parasail_sse41_table OBJECT ${SRC_SSE41} )
    ADD_LIBRARY( parasail_sse41_rowcol OBJECT ${SRC_SSE41} )
    ADD_LIBRARY( parasail_sse41_trace OBJECT ${SRC_TRACE_SSE41} )
    SET_TARGET_PROPERTIES( parasail_sse41 PROPERTIES COMPILE_FLAGS ${SSE41_C_FLAGS} )
    SET_TARGET_PROPERTIES( parasail_sse41_table PROPERTIES COMPILE_FLAGS ${SSE41_C_FLAGS} )
    SET_TARGET_PROPERTIES( parasail_sse41_rowcol PROPERTIES COMPILE_FLAGS ${SSE41_C_FLAGS} )
    SET_TARGET_PROPERTIES( parasail_sse41_trace PROPERTIES COMPILE_FLAGS ${SSE41_C_FLAGS} )
ELSE( )
    ADD_LIBRARY( parasail_sse41 OBJECT cmake/sse41_dummy.c )
    ADD_LIBRARY( parasail_sse41_table OBJECT cmake/sse41_dummy.c )
    ADD_LIBRARY( parasail_sse41_rowcol OBJECT cmake/sse41_dummy.c )
    ADD_LIBRARY( parasail_sse41_trace OBJECT cmake/sse41_dummy.c )
ENDIF( )
SET_TARGET_PROPERTIES( parasail_sse41_table PROPERTIES COMPILE_DEFINITIONS PARASAIL_TABLE )
SET_TARGET_PROPERTIES( parasail_sse41_rowcol PROPERTIES COMPILE_DEFINITIONS PARASAIL_ROWCOL )
SET_TARGET_PROPERTIES( parasail_sse41_trace PROPERTIES COMPILE_DEFINITIONS PARASAIL_TRACE )

IF( AVX2_FOUND )
    ADD_LIBRARY( parasail_avx2 OBJECT ${SRC_AVX2} src/memory_avx2.c )
    ADD_LIBRARY( parasail_avx2_table OBJECT ${SRC_AVX2} )
    ADD_LIBRARY( parasail_avx2_rowcol OBJECT ${SRC_AVX2} )
    ADD_LIBRARY( parasail_avx2_trace OBJECT ${SRC_TRACE_AVX2} )
    SET_TARGET_PROPERTIES( parasail_avx2 PROPERTIES COMPILE_FLAGS ${AVX2_C_FLAGS} )
    SET_TARGET_PROPERTIES( parasail_avx2_table PROPERTIES COMPILE_FLAGS ${AVX2_C_FLAGS} )
    SET_TARGET_PROPERTIES( parasail_avx2_rowcol PROPERTIES COMPILE_FLAGS ${AVX2_C_FLAGS} )
    SET_TARGET_PROPERTIES( parasail_avx2_trace PROPERTIES COMPILE_FLAGS ${AVX2_C_FLAGS} )
ELSE( )
    ADD_LIBRARY( parasail_avx2 OBJECT cmake/avx2_dummy.c )
    ADD_LIBRARY( parasail_avx2_table OBJECT cmake/avx2_dummy.c )
    ADD_LIBRARY( parasail_avx2_rowcol OBJECT cmake/avx2_dummy.c )
    ADD_LIBRARY( parasail_avx2_trace OBJECT cmake/avx2_dummy.c )
ENDIF( )
SET_TARGET_PROPERTIES( parasail_avx2_table PROPERTIES COMPILE_DEFINITIONS PARASAIL_TABLE )
SET_TARGET_PROPERTIES( parasail_avx2_rowcol PROPERTIES COMPILE_DEFINITIONS PARASAIL_ROWCOL )
SET_TARGET_PROPERTIES( parasail_avx2_trace PROPERTIES COMPILE_DEFINITIONS PARASAIL_TRACE )

IF( ALTIVEC_FOUND )
    ADD_LIBRARY( parasail_altivec OBJECT ${SRC_ALTIVEC} src/memory_altivec.c )
    ADD_LIBRARY( parasail_altivec_table OBJECT ${SRC_ALTIVEC} )
    ADD_LIBRARY( parasail_altivec_rowcol OBJECT ${SRC_ALTIVEC} )
    ADD_LIBRARY( parasail_altivec_trace OBJECT ${SRC_TRACE_ALTIVEC} )
    SET_TARGET_PROPERTIES( parasail_altivec PROPERTIES COMPILE_FLAGS ${ALTIVEC_C_FLAGS} )
    SET_TARGET_PROPERTIES( parasail_altivec_table PROPERTIES COMPILE_FLAGS ${ALTIVEC_C_FLAGS} )
    SET_TARGET_PROPERTIES( parasail_altivec_rowcol PROPERTIES COMPILE_FLAGS ${ALTIVEC_C_FLAGS} )
    SET_TARGET_PROPERTIES( parasail_altivec_trace PROPERTIES COMPILE_FLAGS ${ALTIVEC_C_FLAGS} )
ELSE( )
    ADD_LIBRARY( parasail_altivec OBJECT cmake/altivec_dummy.c )
    ADD_LIBRARY( parasail_altivec_table OBJECT cmake/altivec_dummy.c )
    ADD_LIBRARY( parasail_altivec_rowcol OBJECT cmake/altivec_dummy.c )
    ADD_LIBRARY( parasail_altivec_trace OBJECT cmake/altivec_dummy.c )
ENDIF( )
SET_TARGET_PROPERTIES( parasail_altivec_table PROPERTIES COMPILE_DEFINITIONS PARASAIL_TABLE )
SET_TARGET_PROPERTIES( parasail_altivec_rowcol PROPERTIES COMPILE_DEFINITIONS PARASAIL_ROWCOL )
SET_TARGET_PROPERTIES( parasail_altivec_trace PROPERTIES COMPILE_DEFINITIONS PARASAIL_TRACE )

IF( NEON_FOUND )
    ADD_LIBRARY( parasail_neon OBJECT ${SRC_NEON} src/memory_neon.c )
    ADD_LIBRARY( parasail_neon_table OBJECT ${SRC_NEON} )
    ADD_LIBRARY( parasail_neon_rowcol OBJECT ${SRC_NEON} )
    ADD_LIBRARY( parasail_neon_trace OBJECT ${SRC_TRACE_NEON} )
    SET_TARGET_PROPERTIES( parasail_neon PROPERTIES COMPILE_FLAGS ${NEON_C_FLAGS} )
    SET_TARGET_PROPERTIES( parasail_neon_table PROPERTIES COMPILE_FLAGS ${NEON_C_FLAGS} )
    SET_TARGET_PROPERTIES( parasail_neon_rowcol PROPERTIES COMPILE_FLAGS ${NEON_C_FLAGS} )
    SET_TARGET_PROPERTIES( parasail_neon_trace PROPERTIES COMPILE_FLAGS ${NEON_C_FLAGS} )
ELSE( )
    ADD_LIBRARY( parasail_neon OBJECT cmake/neon_dummy.c )
    ADD_LIBRARY( parasail_neon_table OBJECT cmake/neon_dummy.c )
    ADD_LIBRARY( parasail_neon_rowcol OBJECT cmake/neon_dummy.c )
    ADD_LIBRARY( parasail_neon_trace OBJECT cmake/neon_dummy.c )
ENDIF( )
SET_TARGET_PROPERTIES( parasail_neon_table PROPERTIES COMPILE_DEFINITIONS PARASAIL_TABLE )
SET_TARGET_PROPERTIES( parasail_neon_rowcol PROPERTIES COMPILE_DEFINITIONS PARASAIL_ROWCOL )
SET_TARGET_PROPERTIES( parasail_neon_trace PROPERTIES COMPILE_DEFINITIONS PARASAIL_TRACE )

ADD_LIBRARY( sais OBJECT contrib/sais-lite-lcp/sais.c )

#IF( SSE2_FOUND )
#    ADD_LIBRARY( ssw OBJECT contrib/ssw/ssw.c )
#    SET_TARGET_PROPERTIES( ssw PROPERTIES COMPILE_FLAGS ${SSE2_C_FLAGS} )
#ELSE( )
#    ADD_LIBRARY( ssw OBJECT cmake/ssw_dummy.c )
#ENDIF( )

SET( maybe_def )
IF( WIN32 )
    IF( BUILD_SHARED_LIBS )
        LIST( APPEND maybe_def cmake/parasail.def )
    ENDIF( )
    # Disable C4244 (implicit integer truncation) warning on MSVC builds
    # because there are too many caused by the 64-bit wide integer
    # versions of the alignment functions.
    SET( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} /wd4244")
ENDIF( )

ADD_LIBRARY( parasail ${LIB_TYPE}
    $<TARGET_OBJECTS:parasail_core>
    $<TARGET_OBJECTS:parasail_novec>
    $<TARGET_OBJECTS:parasail_novec_table>
    $<TARGET_OBJECTS:parasail_novec_rowcol>
    $<TARGET_OBJECTS:parasail_novec_trace>
    $<TARGET_OBJECTS:parasail_sse2>
    $<TARGET_OBJECTS:parasail_sse2_table>
    $<TARGET_OBJECTS:parasail_sse2_rowcol>
    $<TARGET_OBJECTS:parasail_sse2_trace>
    $<TARGET_OBJECTS:parasail_sse41>
    $<TARGET_OBJECTS:parasail_sse41_table>
    $<TARGET_OBJECTS:parasail_sse41_rowcol>
    $<TARGET_OBJECTS:parasail_sse41_trace>
    $<TARGET_OBJECTS:parasail_avx2>
    $<TARGET_OBJECTS:parasail_avx2_table>
    $<TARGET_OBJECTS:parasail_avx2_rowcol>
    $<TARGET_OBJECTS:parasail_avx2_trace>
    $<TARGET_OBJECTS:parasail_altivec>
    $<TARGET_OBJECTS:parasail_altivec_table>
    $<TARGET_OBJECTS:parasail_altivec_rowcol>
    $<TARGET_OBJECTS:parasail_altivec_trace>
    $<TARGET_OBJECTS:parasail_neon>
    $<TARGET_OBJECTS:parasail_neon_table>
    $<TARGET_OBJECTS:parasail_neon_rowcol>
    $<TARGET_OBJECTS:parasail_neon_trace>
    ${maybe_def}
)

TARGET_INCLUDE_DIRECTORIES(parasail PUBLIC ${CMAKE_CURRENT_LIST_DIR})

IF( ZLIB_FOUND )
    TARGET_LINK_LIBRARIES( parasail ${ZLIB_LIBRARIES} )
    TARGET_INCLUDE_DIRECTORIES( parasail PUBLIC ${ZLIB_INCLUDE_DIR} )
    SET( HAVE_ZLIB TRUE CACHE INTERNAL "" )
ELSE( )
    SET( HAVE_ZLIB FALSE CACHE INTERNAL "" )
ENDIF( )


CHECK_C_SOURCE_COMPILES(
"
int main(void) {
    long long val = 116444736000000000i64;
    return 0;
}
"
    INT64_LITERAL_SUFFIX_I64 )
CHECK_C_SOURCE_COMPILES(
"
int main(void) {
    long long val = 116444736000000000LL;
    return 0;
}
"
    INT64_LITERAL_SUFFIX_LL )

IF( HAVE_CLOCK_GETTIME_NATIVE )
ELSEIF( HAVE_CLOCK_GETTIME_RT )
    TARGET_LINK_LIBRARIES( parasail rt )
ELSEIF( HAVE_CLOCK_GETTIME_PTHREAD )
    TARGET_LINK_LIBRARIES( parasail pthread )
ENDIF( )
# check for monotonic
CHECK_C_SOURCE_COMPILES(
"
#define _POSIX_C_SOURCE 199309L
#include <time.h>
double parasail_time(void) {
    struct timespec ts;
    /* Works on FreeBSD */
    long retval = clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)(ts.tv_sec) + (double)(ts.tv_nsec)/1000000000.0;
}
int main(void) {
    double t = parasail_time();
    return 0;
}
"
    HAVE_CLOCK_GETTIME_MONOTONIC )
# check for realtime
CHECK_C_SOURCE_COMPILES(
"
#define _POSIX_C_SOURCE 199309L
#include <time.h>
double parasail_time(void) {
    struct timespec ts;
    /* Works on FreeBSD */
    long retval = clock_gettime(CLOCK_REALTIME, &ts);
    return (double)(ts.tv_sec) + (double)(ts.tv_nsec)/1000000000.0;
}
int main(void) {
    double t = parasail_time();
    return 0;
}
"
    HAVE_CLOCK_GETTIME_REALTIME )

IF( HAVE_SQRT_NATIVE )
ELSEIF( HAVE_SQRT_M )
    TARGET_LINK_LIBRARIES( parasail m )
ENDIF( )

SET( maybe_getopt )
IF( NOT HAVE_GETOPT AND HAVE_WINDOWS_H )
    LIST( APPEND maybe_getopt contrib/wingetopt/src/getopt.c )
    MESSAGE( STATUS "Using contrib getopt" )
ENDIF( )

ADD_EXECUTABLE( parasail_aligner apps/parasail_aligner.cpp contrib/getMemorySize.c $<TARGET_OBJECTS:sais> ${maybe_getopt} )
IF( OPENMP_FOUND )
    SET_TARGET_PROPERTIES( parasail_aligner PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS} )
    IF( NOT WIN32 )
        SET_TARGET_PROPERTIES( parasail_aligner PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS} )
    ENDIF( )
ENDIF( )
TARGET_LINK_LIBRARIES( parasail_aligner parasail )

ADD_EXECUTABLE( parasail_stats apps/parasail_stats.c )
TARGET_LINK_LIBRARIES( parasail_stats parasail )

ADD_EXECUTABLE( traceback tests/traceback.c ${maybe_getopt} )
TARGET_LINK_LIBRARIES( traceback parasail )

ADD_EXECUTABLE( test_align tests/test_align.c ${maybe_getopt} )
TARGET_LINK_LIBRARIES( test_align parasail )

ADD_EXECUTABLE( test_gcups tests/test_gcups.c ${maybe_getopt} )
IF( OPENMP_FOUND )
    SET_TARGET_PROPERTIES( test_gcups PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS} )
    IF( NOT WIN32 )
        SET_TARGET_PROPERTIES( test_gcups PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS} )
    ENDIF( )
ENDIF( )
TARGET_LINK_LIBRARIES( test_gcups parasail )

ADD_EXECUTABLE( test_isa tests/test_isa.c )
TARGET_LINK_LIBRARIES( test_isa parasail )

ADD_EXECUTABLE( test_matrix tests/test_matrix.c )
TARGET_LINK_LIBRARIES( test_matrix parasail )

INSTALL( FILES parasail.h DESTINATION include )
INSTALL( DIRECTORY parasail DESTINATION include )
INSTALL( TARGETS parasail parasail_aligner parasail_stats
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

SET( CPACK_PACKAGE_VERSION_MAJOR ${MAJOR} )
SET( CPACK_PACKAGE_VERSION_MINOR ${MINOR} )
SET( CPACK_PACKAGE_VERSION_PATCH ${PATCH} )
SET( CPACK_PACKAGE_VERSION ${MAJOR}.${MINOR}.${PATCH} )
SET( CPACK_PACKAGE_NAME "${PROJECT_NAME}" )
SET( CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_LIST_DIR}/README.md" )
SET( CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_LIST_DIR}/LICENSE.md" )
SET( CPACK_OUTPUT_FILE_PREFIX dist )
IF( WIN32 )
    SET( CPACK_GENERATOR "ZIP" )
    IF( "${CMAKE_SIZEOF_VOID_P}" EQUAL "8" )
        SET( CPACK_SYSTEM_NAME "win64" )
    ELSE( )
        SET( CPACK_SYSTEM_NAME "win32" )
    ENDIF( )
    SET( CPACK_PACKAGE_FILE_NAME ${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${CPACK_SYSTEM_NAME}-${CMAKE_VS_PLATFORM_TOOLSET} )
ELSEIF( APPLE )
    SET( CPACK_GENERATOR "TGZ" )
    # there must be an easier way to do this
    EXEC_PROGRAM( "/usr/bin/sw_vers" OUTPUT_VARIABLE SW_VERS )
    STRING( REGEX REPLACE "BuildVersion:.*" "" SW_VERS ${SW_VERS} )
    STRING( REGEX REPLACE "[ \t\n]+" "" SW_VERS ${SW_VERS} )
    STRING( REPLACE "ProductName:MacOSXProductVersion:" "" SW_VERS ${SW_VERS} )
    SET( CPACK_PACKAGE_FILE_NAME ${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${CMAKE_SYSTEM_NAME}-${SW_VERS} )
ELSE( )
    SET( CPACK_GENERATOR "TGZ" )
ENDIF( )

INCLUDE( CPack )
