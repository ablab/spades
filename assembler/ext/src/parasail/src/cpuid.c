/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 *
 * These functions borrow heavily, if not verbatim, from the files in
 * the parsail contrib directory. Please see contrib for details and
 * copyrights of the respective authors/sources.
 */
#include "config.h"

#include "parasail/parasail/cpuid.h"

#include <stdint.h>
#if defined(_MSC_VER)
# include <intrin.h>
#endif

static void run_cpuid(uint32_t eax, uint32_t ecx, uint32_t* abcd)
{
#if defined(_MSC_VER)
    __cpuidex(abcd, eax, ecx);
#else
    uint32_t ebx=0, edx=0;
# if defined( __i386__ ) && defined ( __PIC__ )
     /* in case of PIC under 32-bit EBX cannot be clobbered */
    __asm__ ( "movl %%ebx, %%edi \n\t cpuid \n\t xchgl %%ebx, %%edi" : "=D" (ebx),
# else
    __asm__ ( "cpuid" : "+b" (ebx),
# endif
              "+a" (eax), "+c" (ecx), "=d" (edx) );
    abcd[0] = eax; abcd[1] = ebx; abcd[2] = ecx; abcd[3] = edx;
#endif
}     


#if defined(__INTEL_COMPILER) && (__INTEL_COMPILER >= 1300)

#include <immintrin.h>

static int check_4th_gen_intel_core_features()
{
    const int the_4th_gen_features = 
        (_FEATURE_AVX2 | _FEATURE_FMA | _FEATURE_BMI | _FEATURE_LZCNT | _FEATURE_MOVBE);
    return _may_i_use_cpu_feature( the_4th_gen_features );
}

static int has_intel_avx512f_features()
{
    const unsigned long avx512f_features =
        (_FEATURE_AVX512F | _FEATURE_AVX512ER | _FEATURE_AVX512PF | _FEATURE_AVX512CD);
    return _may_i_use_cpu_feature( avx512f_features );
}

static int has_intel_avx512bw_features()
{
    const unsigned long avx512f_features =
        (_FEATURE_AVX512VL | _FEATURE_AVX512BW | _FEATURE_AVX512DQ);
    return _may_i_use_cpu_feature( avx512f_features );
}

static int has_intel_avx512vbmi_features()
{
    const unsigned long avx512vbmi_features = (_FEATURE_AVX512VBMI);
    return _may_i_use_cpu_feature( avx512vbmi_features );
}

#else /* non-Intel compiler */

static int check_xcr0_ymm() 
{
#if HAVE_XGETBV
    uint32_t xcr0;
#if defined(_MSC_VER)
    xcr0 = (uint32_t)_xgetbv(0);  /* min VS2010 SP1 compiler is required */
#else
    __asm__ ("xgetbv" : "=a" (xcr0) : "c" (0) : "%edx" );
#endif
    return ((xcr0 & 6) == 6); /* checking if xmm and ymm state are enabled in XCR0 */
#else
    return 0;
#endif
}

static int check_xcr0_zmm()
{
#if HAVE_XGETBV
    uint32_t xcr0;
    uint32_t zmm_ymm_xmm = (7U << 5) | (1U << 2) | (1U << 1);
#if defined(_MSC_VER)
    xcr0 = (uint32_t)_xgetbv(0);  /* min VS2010 SP1 compiler is required */
#else
    __asm__ ("xgetbv" : "=a" (xcr0) : "c" (0) : "%edx" );
#endif
    return ((xcr0 & zmm_ymm_xmm) == zmm_ymm_xmm); /* check if xmm, zmm and zmm state are enabled in XCR0 */
#else
    return 0;
#endif
}

static int check_4th_gen_intel_core_features()
{
    uint32_t abcd[4];
    uint32_t fma_movbe_osxsave_mask = ((1U << 12) | (1U << 22) | (1U << 27));
    uint32_t avx2_bmi12_mask = (1U << 5) | (1U << 3) | (1U << 8);

    /* CPUID.(EAX=01H, ECX=0H):ECX.FMA[bit 12]==1   && 
       CPUID.(EAX=01H, ECX=0H):ECX.MOVBE[bit 22]==1 && 
       CPUID.(EAX=01H, ECX=0H):ECX.OSXSAVE[bit 27]==1 */
    run_cpuid( 1, 0, abcd );
    if ( (abcd[2] & fma_movbe_osxsave_mask) != fma_movbe_osxsave_mask ) 
        return 0;

    if ( ! check_xcr0_ymm() )
        return 0;

    /*  CPUID.(EAX=07H, ECX=0H):EBX.AVX2[bit 5]==1  &&
        CPUID.(EAX=07H, ECX=0H):EBX.BMI1[bit 3]==1  &&
        CPUID.(EAX=07H, ECX=0H):EBX.BMI2[bit 8]==1  */
    run_cpuid( 7, 0, abcd );
    if ( (abcd[1] & avx2_bmi12_mask) != avx2_bmi12_mask ) 
        return 0;

    /* CPUID.(EAX=80000001H):ECX.LZCNT[bit 5]==1 */
    run_cpuid( 0x80000001, 0, abcd );
    if ( (abcd[2] & (1U << 5)) == 0)
        return 0;

    return 1;
}

static int has_intel_avx512f_features() {
    uint32_t abcd[4];
    uint32_t osxsave_mask = (1U << 27); /* OSX. */
    uint32_t avx512f_mask = (1U << 16) | /* AVX-512F */
        (1U << 26) | /* AVX-512PF */
        (1U << 27) | /* AVX-512ER */
        (1U << 28);  /* AVX-512CD */

    /* CPUID.(EAX=01H, ECX=0H):ECX.OSXSAVE[bit 27]==1 */
    run_cpuid( 1, 0, abcd );
    if ( (abcd[2] & osxsave_mask) != osxsave_mask ) 
        return 0;

    if ( ! check_xcr0_zmm() )
        return 0;

    /*  CPUID.(EAX=07H, ECX=0H):EBX.AVX-512F [bit 16]==1  &&
        CPUID.(EAX=07H, ECX=0H):EBX.AVX-512PF[bit 26]==1  &&
        CPUID.(EAX=07H, ECX=0H):EBX.AVX-512ER[bit 27]==1  &&
        CPUID.(EAX=07H, ECX=0H):EBX.AVX-512CD[bit 28]==1  */
    run_cpuid( 7, 0, abcd );
    if ( (abcd[1] & avx512f_mask) != avx512f_mask ) 
        return 0;

    return 1;
}

static int has_intel_avx512bw_features() {
    uint32_t abcd[4];
    uint32_t osxsave_mask = (1U << 27); /* OSX. */
    uint32_t avx512bw_mask = (1U << 30) | /* AVX-512BW */
        (1U << 17) | /* AVX-512DQ */
        (1U << 31);  /* AVX-512VL */

    /* CPUID.(EAX=01H, ECX=0H):ECX.OSXSAVE[bit 27]==1 */
    run_cpuid( 1, 0, abcd );
    if ( (abcd[2] & osxsave_mask) != osxsave_mask ) 
        return 0;

    if ( ! check_xcr0_zmm() )
        return 0;

    /*  CPUID.(EAX=07H, ECX=0H):EBX.AVX-512BW[bit 30]==1  &&
        CPUID.(EAX=07H, ECX=0H):EBX.AVX-512DQ[bit 17]==1  &&
        CPUID.(EAX=07H, ECX=0H):EBX.AVX-512VL[bit 31]==1  */
    run_cpuid( 7, 0, abcd );
    if ( (abcd[1] & avx512bw_mask) != avx512bw_mask ) 
        return 0;

    return 1;
}

static int has_intel_avx512vbmi_features()
{
    uint32_t abcd[4];
    uint32_t osxsave_mask = (1U << 27); /* OSX. */
    uint32_t avx512vbmi_mask = (1U << 1); /* AVX-512VBMI */

    /* CPUID.(EAX=01H, ECX=0H):ECX.OSXSAVE[bit 27]==1 */
    run_cpuid( 1, 0, abcd );
    if ( (abcd[2] & osxsave_mask) != osxsave_mask ) 
        return 0;

    if ( ! check_xcr0_zmm() )
        return 0;

    /*  CPUID.(EAX=07H, ECX=0H):ECX.AVX-512VBMI[bit 1]==1 */
    run_cpuid( 7, 0, abcd );
    if ( (abcd[2] & avx512vbmi_mask) != avx512vbmi_mask ) 
        return 0;

    return 1;
}

#endif /* non-Intel compiler */


static int check_sse41()
{
    uint32_t info[4];
    uint32_t nIds;

    run_cpuid(0, 0, info);
    nIds = info[0];

    /*  Detect Instruction Set */
    if (nIds >= 1){
        run_cpuid(0x00000001, 0, info);
        return (info[2] & (1U << 19)) != 0;
    }

    return 0;
}


static int check_sse2()
{
    uint32_t info[4];
    uint32_t nIds;

    run_cpuid(0, 0, info);
    nIds = info[0];

    /*  Detect Instruction Set */
    if (nIds >= 1){
        run_cpuid(0x00000001, 0, info);
        return (info[3] & (1U << 26)) != 0;
    }

    return 0;
}


int parasail_can_use_avx512vbmi()
{
    static int avx512vbmi_features_available = -1;
    /* test is performed once */
    if (avx512vbmi_features_available < 0 )
        avx512vbmi_features_available = has_intel_avx512vbmi_features();
    return avx512vbmi_features_available;
}

int parasail_can_use_avx512bw()
{
    static int avx512bw_features_available = -1;
    /* test is performed once */
    if (avx512bw_features_available < 0 )
        avx512bw_features_available = has_intel_avx512bw_features();
    return avx512bw_features_available;
}

int parasail_can_use_avx512f()
{
    static int avx512f_features_available = -1;
    /* test is performed once */
    if (avx512f_features_available < 0 )
        avx512f_features_available = has_intel_avx512f_features();
    return avx512f_features_available;
}

int parasail_can_use_avx2()
{
    static int the_4th_gen_features_available = -1;
    /* test is performed once */
    if (the_4th_gen_features_available < 0 )
        the_4th_gen_features_available = check_4th_gen_intel_core_features();

    return the_4th_gen_features_available;
}

int parasail_can_use_sse41()
{
    static int can_use_sse41 = -1;
    /* test is performed once */
    if (can_use_sse41 < 0)
        can_use_sse41 = check_sse41();

    return can_use_sse41;
}

int parasail_can_use_sse2()
{
    static int can_use_sse2 = -1;
    /* test is performed once */
    if (can_use_sse2 < 0)
        can_use_sse2 = check_sse2();

    return can_use_sse2;
}

int parasail_can_use_altivec()
{
    return 0;
}

int parasail_can_use_neon()
{
    return 0;
}

