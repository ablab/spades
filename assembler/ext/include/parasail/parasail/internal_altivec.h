/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyb (c) 2015 Battelle Memorial Institute.
 */
#ifndef _PARASAIL_INTERNAL_ALTIVEC_H_
#define _PARASAIL_INTERNAL_ALTIVEC_H_

#include <stdint.h>

#include <altivec.h>

#ifdef __cplusplus
extern "C" {
#endif

#if defined(_MSC_VER)
#define ALIGNED16 __declspec( align( 16 ) )
#else
#define ALIGNED16 __attribute__ ((__aligned__ (16)))
#endif

typedef ALIGNED16 vector signed char      vec128i;
typedef ALIGNED16 vector unsigned char    vec16ub;
typedef ALIGNED16 vector signed char      vec16sb;
typedef ALIGNED16 vector signed short     vec8sh;
typedef ALIGNED16 vector signed int       vec4sw;
typedef ALIGNED16 vector signed long long vec2sd;

typedef ALIGNED16 union vec128i_8 {
    vec128i m;
    int8_t v[16];
} vec128i_8_t;

typedef ALIGNED16 union vec128i_16 {
    vec128i m;
    int16_t v[8];
} vec128i_16_t;

typedef ALIGNED16 union vec128i_32 {
    vec128i m;
    int32_t v[4];
} vec128i_32_t;

typedef ALIGNED16 union vec128i_64 {
    vec128i m;
    int64_t v[2];
} vec128i_64_t;

extern vec128i * parasail_memalign_vec128i(size_t alignment, size_t size);

extern void parasail_memset_vec128i(vec128i *b, vec128i c, size_t len);

extern void parasail_free_vec128i(void *ptr);

static inline vec128i _mm_add_epi8(vec128i a, vec128i b)
{
    return (vec128i) vec_add((vec16sb) a, (vec16sb) b);
}

static inline vec128i _mm_add_epi16(vec128i a, vec128i b)
{
    return (vec128i) vec_add((vec8sh) a, (vec8sh) b);
}

static inline vec128i _mm_add_epi32(vec128i a, vec128i b)
{
    return (vec128i) vec_add((vec4sw) a, (vec4sw) b);
}

static inline vec128i _mm_add_epi64(vec128i a, vec128i b)
{
    return (vec128i) vec_add((vec2sd) a, (vec2sd) b);
}

static inline vec128i _mm_adds_epi8(vec128i a, vec128i b)
{
    return (vec128i) vec_adds((vec16sb) a, (vec16sb) b);
}

static inline vec128i _mm_adds_epi16(vec128i a, vec128i b)
{
    return (vec128i) vec_adds((vec8sh) a, (vec8sh) b);
}

static inline vec128i _mm_sub_epi8(vec128i a, vec128i b)
{
    return (vec128i) vec_sub((vec16sb) a, (vec16sb) b);
}

static inline vec128i _mm_sub_epi16(vec128i a, vec128i b)
{
    return (vec128i) vec_sub((vec8sh) a, (vec8sh) b);
}

static inline vec128i _mm_sub_epi32(vec128i a, vec128i b)
{
    return (vec128i) vec_sub((vec4sw) a, (vec4sw) b);
}

static inline vec128i _mm_sub_epi64(vec128i a, vec128i b)
{
    return (vec128i) vec_sub((vec2sd) a, (vec2sd) b);
}

static inline vec128i _mm_subs_epi8(vec128i a, vec128i b)
{
    return (vec128i) vec_subs((vec16sb) a, (vec16sb) b);
}

static inline vec128i _mm_subs_epi16(vec128i a, vec128i b)
{
    return (vec128i) vec_subs((vec8sh) a, (vec8sh) b);
}

static inline vec128i _mm_cmpeq_epi8(vec128i a, vec128i b)
{
    return (vec128i) vec_cmpeq((vec16sb) a, (vec16sb) b);
}

static inline vec128i _mm_cmpeq_epi16(vec128i a, vec128i b)
{
    return (vec128i) vec_cmpeq((vec8sh) a, (vec8sh) b);
}

static inline vec128i _mm_cmpeq_epi32(vec128i a, vec128i b)
{
    return (vec128i) vec_cmpeq((vec4sw) a, (vec4sw) b);
}

static inline vec128i _mm_cmpeq_epi64(vec128i a, vec128i b)
{
    return (vec128i) vec_cmpeq((vec2sd) a, (vec2sd) b);
}

static inline vec128i _mm_cmpgt_epi8(vec128i a, vec128i b)
{
    return (vec128i) vec_cmpgt((vec16sb) a, (vec16sb) b);
}

static inline vec128i _mm_cmpgt_epi16(vec128i a, vec128i b)
{
    return (vec128i) vec_cmpgt((vec8sh) a, (vec8sh) b);
}

static inline vec128i _mm_cmpgt_epi32(vec128i a, vec128i b)
{
    return (vec128i) vec_cmpgt((vec4sw) a, (vec4sw) b);
}

static inline vec128i _mm_cmpgt_epi64(vec128i a, vec128i b)
{
    return (vec128i) vec_cmpgt((vec2sd) a, (vec2sd) b);
}

static inline vec128i _mm_cmplt_epi8(vec128i a, vec128i b)
{
    return (vec128i) vec_cmplt((vec16sb) a, (vec16sb) b);
}

static inline vec128i _mm_cmplt_epi16(vec128i a, vec128i b)
{
    return (vec128i) vec_cmplt((vec8sh) a, (vec8sh) b);
}

static inline vec128i _mm_cmplt_epi32(vec128i a, vec128i b)
{
    return (vec128i) vec_cmplt((vec4sw) a, (vec4sw) b);
}

static inline vec128i _mm_cmplt_epi64(vec128i a, vec128i b)
{
    return (vec128i) vec_cmplt((vec2sd) a, (vec2sd) b);
}

static inline int8_t _mm_extract_epi8(vec128i v, int imm)
{
    vec128i_8_t A;
    A.m = v;
#ifdef WORDS_BIGENDIAN
    return A.v[15-imm];
#else
    return A.v[imm];
#endif
}

static inline int16_t _mm_extract_epi16(vec128i v, int imm)
{
    vec128i_16_t A;
    A.m = v;
#ifdef WORDS_BIGENDIAN
    return A.v[7-imm];
#else
    return A.v[imm];
#endif
}

static inline int32_t _mm_extract_epi32(vec128i v, int imm)
{
    vec128i_32_t A;
    A.m = v;
#ifdef WORDS_BIGENDIAN
    return A.v[3-imm];
#else
    return A.v[imm];
#endif
}

static inline int64_t _mm_extract_epi64(vec128i v, int imm)
{
    vec128i_64_t A;
    A.m = v;
#ifdef WORDS_BIGENDIAN
    return A.v[1-imm];
#else
    return A.v[imm];
#endif
}

static inline vec128i _mm_insert_epi8(vec128i v, int8_t value, int imm)
{
    vec128i_8_t A;
    A.m = v;
#ifdef WORDS_BIGENDIAN
    A.v[15-imm] = value;
#else
    A.v[imm] = value;
#endif
    return A.m;
}

static inline vec128i _mm_insert_epi16(vec128i v, int16_t value, int imm)
{
    vec128i_16_t A;
    A.m = v;
#ifdef WORDS_BIGENDIAN
    A.v[7-imm] = value;
#else
    A.v[imm] = value;
#endif
    return A.m;
}

static inline vec128i _mm_insert_epi32(vec128i v, int32_t value, int imm)
{
    vec128i_32_t A;
    A.m = v;
#ifdef WORDS_BIGENDIAN
    A.v[3-imm] = value;
#else
    A.v[imm] = value;
#endif
    return A.m;
}

static inline vec128i _mm_insert_epi64(vec128i v, int64_t value, int imm)
{
    vec128i_64_t A;
    A.m = v;
#ifdef WORDS_BIGENDIAN
    A.v[1-imm] = value;
#else
    A.v[imm] = value;
#endif
    return A.m;
}

static inline vec128i _mm_load_si128(const vec128i *address)
{
    return (vec128i) vec_ld(0, (vec16ub*) address);
}

static inline void _mm_store_si128(vec128i *address, vec128i v)
{
    vec_st(v, 0, address);
}

static inline vec128i _mm_max_epi8(vec128i a, vec128i b)
{
    return (vec128i) vec_max((vec16sb) a, (vec16sb) b);
}

static inline vec128i _mm_max_epi16(vec128i a, vec128i b)
{
    return (vec128i) vec_max((vec8sh) a, (vec8sh) b);
}

static inline vec128i _mm_max_epi32(vec128i a, vec128i b)
{
    return (vec128i) vec_max((vec4sw) a, (vec4sw) b);
}

static inline vec128i _mm_max_epi64(vec128i a, vec128i b)
{
    return (vec128i) vec_max((vec2sd) a, (vec2sd) b);
}

static inline vec128i _mm_min_epi8(vec128i a, vec128i b)
{
    return (vec128i) vec_min((vec16sb) a, (vec16sb) b);
}

static inline vec128i _mm_min_epi16(vec128i a, vec128i b)
{
    return (vec128i) vec_min((vec8sh) a, (vec8sh) b);
}

static inline vec128i _mm_min_epi32(vec128i a, vec128i b)
{
    return (vec128i) vec_min((vec4sw) a, (vec4sw) b);
}

static inline vec128i _mm_min_epi64(vec128i a, vec128i b)
{
    return (vec128i) vec_min((vec2sd) a, (vec2sd) b);
}

static inline int _mm_movemask_epi8(vec128i v)
{
    int ret = 0;
    vec128i_8_t a;
    a.m = v;
#ifdef WORDS_BIGENDIAN
    ret |= (a.v[0]  & 0x80) << (15-7);
    ret |= (a.v[1]  & 0x80) << (14-7);
    ret |= (a.v[2]  & 0x80) << (13-7);
    ret |= (a.v[3]  & 0x80) << (12-7);
    ret |= (a.v[4]  & 0x80) << (11-7);
    ret |= (a.v[5]  & 0x80) << (10-7);
    ret |= (a.v[6]  & 0x80) <<  (9-7);
    ret |= (a.v[7]  & 0x80) <<  (8-7);
    ret |= (a.v[8]  & 0x80);
    ret |= (a.v[9]  & 0x80) >>  (7-6);
    ret |= (a.v[10] & 0x80) >>  (7-5);
    ret |= (a.v[11] & 0x80) >>  (7-4);
    ret |= (a.v[12] & 0x80) >>  (7-3);
    ret |= (a.v[13] & 0x80) >>  (7-2);
    ret |= (a.v[14] & 0x80) >>  (7-1);
    ret |= (a.v[15] & 0x80) >>   7;
#else
    ret |= (a.v[15] & 0x80) << (15-7);
    ret |= (a.v[14] & 0x80) << (14-7);
    ret |= (a.v[13] & 0x80) << (13-7);
    ret |= (a.v[12] & 0x80) << (12-7);
    ret |= (a.v[11] & 0x80) << (11-7);
    ret |= (a.v[10] & 0x80) << (10-7);
    ret |= (a.v[9]  & 0x80) <<  (9-7);
    ret |= (a.v[8]  & 0x80) <<  (8-7);
    ret |= (a.v[7]  & 0x80);
    ret |= (a.v[6]  & 0x80) >>  (7-6);
    ret |= (a.v[5]  & 0x80) >>  (7-5);
    ret |= (a.v[4]  & 0x80) >>  (7-4);
    ret |= (a.v[3]  & 0x80) >>  (7-3);
    ret |= (a.v[2]  & 0x80) >>  (7-2);
    ret |= (a.v[1]  & 0x80) >>  (7-1);
    ret |= (a.v[0]  & 0x80) >>   7;
#endif
    return ret;
}

static inline vec128i _mm_packs_epi16(vec128i a, vec128i b)
{
#ifdef WORDS_BIGENDIAN
    return (vec128i) vec_packs((vec8sh) b, (vec8sh) a);
#else
    return (vec128i) vec_packs((vec8sh) a, (vec8sh) b);
#endif
}

static inline vec128i _mm_unpacklo_epi8(vec128i a, vec128i b)
{
    static const vec16ub permute_selector = {
#ifdef WORDS_BIGENDIAN
        0x18, 0x08, 0x19, 0x09, 0x1A, 0x0A, 0x1B, 0x0B, 0x1C, 0x0C, 0x1D, 0x0D, 0x1E, 0x0E, 0x1F, 0x0F
#else
        0x00, 0x10, 0x01, 0x11, 0x02, 0x12, 0x03, 0x13, 0x04, 0x14, 0x05, 0x15, 0x06, 0x16, 0x07, 0x17
#endif
    };
    return vec_perm(a, b, permute_selector);
}

static inline vec128i _mm_unpackhi_epi8(vec128i a, vec128i b)
{
    static const vec16ub permute_selector = {
#ifdef WORDS_BIGENDIAN
        0x10, 0x00, 0x11, 0x01, 0x12, 0x02, 0x13, 0x03, 0x14, 0x04, 0x15, 0x05, 0x16, 0x06, 0x17, 0x07
#else
        0x08, 0x18, 0x09, 0x19, 0x0A, 0x1A, 0x0B, 0x1B, 0x0C, 0x1C, 0x0D, 0x1D, 0x0E, 0x1E, 0x0F, 0x1F
#endif
    };
    return vec_perm(a, b, permute_selector);

}

static inline vec128i _mm_and_si128(vec128i a, vec128i b)
{
    return (vec128i) vec_and((vec16ub) a, (vec16ub) b);
}

static inline vec128i _mm_andnot_si128(vec128i a, vec128i b)
{
    return (vec128i) vec_andc((vec16ub) b, (vec16ub) a);
}

static inline vec128i _mm_or_si128(vec128i a, vec128i b)
{
    return (vec128i) vec_or((vec16ub) a, (vec16ub) b);
}

static inline vec128i _mm_xor_si128(vec128i a, vec128i b)
{
    return (vec128i) vec_xor((vec16ub) a, (vec16ub) b);
}

static inline vec128i _mm_slli_si128(vec128i v, int imm)
{
    if ((unsigned long) imm >= 16) {
        /* SSE2 shifts >= element_size or < 0 produce 0;
         * Altivec/MMX shifts by imm%element_size. */
        return (vec128i) vec_splats (0);
    } else if (imm == 0) {
        return v;
    } else {
        /* The PowerPC byte shift count must be multiplied by 8. */
        /* It need not but can be replicated, which handles both LE and
         * BE shift count positioning. */
        vec128i replicated_count;
        replicated_count = (vec128i) vec_splats ((char)(imm << 3));
        /* AT gcc v7.1 may miscompile vec_sro as vec_slo? */
        return (vec128i) vec_slo (v, replicated_count);
    }
}

static inline vec128i _mm_srli_si128(vec128i v, int imm)
{
    if ((unsigned long) imm >= 16) {
        /* SSE2 shifts >= element_size or < 0 produce 0;
         * Altivec/MMX shifts by imm%element_size. */
        return (vec128i) vec_splats (0);
    } else if (imm == 0) {
        return v;
    } else {
        /* The PowerPC byte shift count must be multiplied by 8. */
        /* It need not but can be replicated, which handles both LE and
         * BE shift count positioning. */
        vec128i replicated_count;
        replicated_count = (vec128i) vec_splats ((char)(imm << 3));
        /* AT gcc v7.1 may miscompile vec_sro as vec_slo? */
        return (vec128i) vec_sro (v, replicated_count);
    }
}

static inline vec128i _mm_rlli_si128(vec128i v, int imm)
{
#if WORDS_BIGENDIAN
    return _mm_or_si128(
            _mm_srli_si128(v,imm),
            _mm_slli_si128(v,16-imm));
#else
    return _mm_or_si128(
            _mm_slli_si128(v,imm),
            _mm_srli_si128(v,16-imm));
#endif
}

static inline vec128i _mm_setzero_si128()
{
    return (vec128i) vec_splats((signed char)0);
}

static inline vec128i _mm_set1_epi8(int8_t scalar)
{
    return (vec128i) vec_splats((signed char)scalar);
}

static inline vec128i _mm_set1_epi16(int16_t scalar)
{
    return (vec128i) vec_splats((signed short)scalar);
}

static inline vec128i _mm_set1_epi32(int32_t scalar)
{
    return (vec128i) vec_splats((signed int)scalar);
}

static inline vec128i _mm_set1_epi64(int64_t scalar)
{
    return (vec128i) vec_splats((signed long long)scalar);
}

static inline vec128i _mm_set_epi8(int8_t c15, int8_t c14, int8_t c13, int8_t c12, int8_t c11, int8_t c10, int8_t c9, int8_t c8, int8_t c7, int8_t c6, int8_t c5, int8_t c4, int8_t c3, int8_t c2, int8_t c1, int8_t c0)
{
    vec128i_8_t a;
#ifdef WORDS_BIGENDIAN
    a.v[0]  = c15;
    a.v[1]  = c14;
    a.v[2]  = c13;
    a.v[3]  = c12;
    a.v[4]  = c11;
    a.v[5]  = c10;
    a.v[6]  = c9;
    a.v[7]  = c8;
    a.v[8]  = c7;
    a.v[9]  = c6;
    a.v[10] = c5;
    a.v[11] = c4;
    a.v[12] = c3;
    a.v[13] = c2;
    a.v[14] = c1;
    a.v[15] = c0;
#else
    a.v[0]  = c0;
    a.v[1]  = c1;
    a.v[2]  = c2;
    a.v[3]  = c3;
    a.v[4]  = c4;
    a.v[5]  = c5;
    a.v[6]  = c6;
    a.v[7]  = c7;
    a.v[8]  = c8;
    a.v[9]  = c9;
    a.v[10] = c10;
    a.v[11] = c11;
    a.v[12] = c12;
    a.v[13] = c13;
    a.v[14] = c14;
    a.v[15] = c15;
#endif
    return (vec128i) a.m;
}

static inline vec128i _mm_set_epi16(int16_t s7, int16_t s6, int16_t s5, int16_t s4, int16_t s3, int16_t s2, int16_t s1, int16_t s0)
{
    vec128i_16_t a;
#if WORDS_BIGENDIAN
    a.v[0] = s7;
    a.v[1] = s6;
    a.v[2] = s5;
    a.v[3] = s4;
    a.v[4] = s3;
    a.v[5] = s2;
    a.v[6] = s1;
    a.v[7] = s0;
#else
    a.v[0] = s0;
    a.v[1] = s1;
    a.v[2] = s2;
    a.v[3] = s3;
    a.v[4] = s4;
    a.v[5] = s5;
    a.v[6] = s6;
    a.v[7] = s7;
#endif
    return (vec128i) a.m;
}

static inline vec128i _mm_set_epi32(int32_t i3, int32_t i2, int32_t i1, int32_t i0)
{
    vec128i_32_t a;
#ifdef WORDS_BIGENDIAN
    a.v[0] = i3;
    a.v[1] = i2;
    a.v[2] = i1;
    a.v[3] = i0;
#else
    a.v[0] = i0;
    a.v[1] = i1;
    a.v[2] = i2;
    a.v[3] = i3;
#endif
    return (vec128i) a.m;
}

static inline vec128i _mm_set_epi64(int64_t l1, int64_t l0)
{
    vec128i_64_t a;
#ifdef WORDS_BIGENDIAN
    a.v[0] = l1;
    a.v[1] = l0;
#else
    a.v[0] = l0;
    a.v[1] = l1;
#endif
    return (vec128i) a.m;
}

static inline int8_t _mm_hmax_epi8(vec128i v)
{
    v = _mm_max_epi8(v, _mm_srli_si128(v, 8));
    v = _mm_max_epi8(v, _mm_srli_si128(v, 4));
    v = _mm_max_epi8(v, _mm_srli_si128(v, 2));
    v = _mm_max_epi8(v, _mm_srli_si128(v, 1));
    return _mm_extract_epi8(v, 0);
}

static inline int16_t _mm_hmax_epi16(vec128i v)
{
    v = _mm_max_epi16(v, _mm_srli_si128(v, 8));
    v = _mm_max_epi16(v, _mm_srli_si128(v, 4));
    v = _mm_max_epi16(v, _mm_srli_si128(v, 2));
    return _mm_extract_epi16(v, 0);
}

static inline int32_t _mm_hmax_epi32(vec128i v)
{
    v = _mm_max_epi32(v, _mm_srli_si128(v, 8));
    v = _mm_max_epi32(v, _mm_srli_si128(v, 4));
    return _mm_extract_epi32(v, 0);
}

static inline int64_t _mm_hmax_epi64(vec128i v)
{
    v = _mm_max_epi64(v, _mm_srli_si128(v, 8));
    return _mm_extract_epi64(v, 0);
}

static inline vec128i _mm_blendv_epi8(vec128i a, vec128i b, vec128i mask)
{
    a = _mm_andnot_si128(mask, a);
    a = _mm_or_si128(a, _mm_and_si128(mask, b));
    return a;
}

#ifdef __cplusplus
}
#endif

#endif /* _PARASAIL_INTERNAL_ALTIVEC_H_ */
