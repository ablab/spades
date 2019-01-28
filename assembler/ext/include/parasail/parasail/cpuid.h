/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#ifndef _PARASAIL_CPUID_H_
#define _PARASAIL_CPUID_H_

#include "parasail/parasail.h"

#ifdef __cplusplus
extern "C" {
#endif

extern int parasail_can_use_avx512vbmi();
extern int parasail_can_use_avx512bw();
extern int parasail_can_use_avx512f();
extern int parasail_can_use_avx2();
extern int parasail_can_use_sse41();
extern int parasail_can_use_sse2();
extern int parasail_can_use_altivec();
extern int parasail_can_use_neon();

#ifdef __cplusplus
}
#endif

#endif /* _PARASAIL_CPUID_H_ */

