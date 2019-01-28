#include "config.h"

#include <stdio.h>
#include <stdlib.h>

#include "parasail/parasail/internal_altivec.h"

static inline void print_8(vec128i v)
{
    vec128i_8_t a;
    a.m = v;
    printf( "%d" ",%d" ",%d" ",%d"
           ",%d" ",%d" ",%d" ",%d"
           ",%d" ",%d" ",%d" ",%d"
           ",%d" ",%d" ",%d" ",%d",
           a.v[0], a.v[1], a.v[2], a.v[3],
           a.v[4], a.v[5], a.v[6], a.v[7],
           a.v[8], a.v[9], a.v[10], a.v[11],
           a.v[12], a.v[13], a.v[14], a.v[15] );
}

static inline void print_16(vec128i v)
{
    vec128i_16_t a;
    a.m = v;
    printf( "%d" ",%d" ",%d" ",%d"
           ",%d" ",%d" ",%d" ",%d",
           a.v[0], a.v[1], a.v[2], a.v[3],
           a.v[4], a.v[5], a.v[6], a.v[7]);
}

static inline void print_32(vec128i v)
{
    vec128i_32_t a;
    a.m = v;
    printf( "%d" ",%d" ",%d" ",%d",
           a.v[0], a.v[1], a.v[2], a.v[3]);
}

static inline void print_64(vec128i v)
{
    vec128i_64_t a;
    a.m = v;
    printf( "%lld" ",%lld",
           a.v[0], a.v[1]);
}

int main(int argc, char **argv)
{
    vec128i vZero8;
    vec128i vZero16;
    vec128i vZero32;
    vec128i vZero64;
    vec128i vOne8;
    vec128i vOne16;
    vec128i vOne32;
    vec128i vOne64;
    vec128i vResult8;
    vec128i vResult16;
    vec128i vResult32;
    vec128i vResult64;
    vec128i vMask8;
    vec128i vMask16;
    vec128i vMask32;
    vec128i vMask64;
    vec128i vEnum8;
    vec128i vEnum16;
    vec128i vEnum32;
    vec128i vEnum64;
    vec128i vEnum8b;
    vec128i vEnum16b;
    vec128i vEnum32b;
    vec128i vEnum64b;
    vec128i *vPointer = NULL;

    vZero8 = _mm_setzero_si128();
    vZero16 = _mm_setzero_si128();
    vZero32 = _mm_setzero_si128();
    vZero64 = _mm_setzero_si128();
    vOne8 = _mm_set1_epi8(1);
    vOne16 = _mm_set1_epi16(1);
    vOne32 = _mm_set1_epi32(1);
    vOne64 = _mm_set1_epi64(1);
    vEnum8 = _mm_set_epi8(15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0);
    vEnum16 = _mm_set_epi16(7,6,5,4,3,2,1,0);
    vEnum32 = _mm_set_epi32(3,2,1,0);
    vEnum64 = _mm_set_epi64(1,0);
    vEnum8b = _mm_set_epi8(-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0);
    vEnum16b = _mm_set_epi16(-7,-6,-5,-4,-3,-2,-1,0);
    vEnum32b = _mm_set_epi32(-3,-2,-1,-0);
    vEnum64b = _mm_set_epi64(-1,-0);
    vPointer = malloc(sizeof(vec128i)*4);

    printf("vZero=");print_8(vZero8);printf("\n");
    printf("vZero=");print_16(vZero16);printf("\n");
    printf("vZero=");print_32(vZero32);printf("\n");
    printf("vZero=");print_64(vZero64);printf("\n");

    printf("vOne=");print_8(vOne8);printf("\n");
    printf("vOne=");print_16(vOne16);printf("\n");
    printf("vOne=");print_32(vOne32);printf("\n");
    printf("vOne=");print_64(vOne64);printf("\n");

    printf("vEnum=");print_8(vEnum8);printf("\n");
    printf("vEnum=");print_16(vEnum16);printf("\n");
    printf("vEnum=");print_32(vEnum32);printf("\n");
    printf("vEnum=");print_64(vEnum64);printf("\n");

    vResult8 = _mm_add_epi8(vOne8, vOne8);
    vResult16 = _mm_add_epi16(vOne16, vOne16);
    vResult32 = _mm_add_epi32(vOne32, vOne32);
    vResult64 = _mm_add_epi64(vOne64, vOne64);

    printf("add=");print_8(vResult8);printf("\n");
    printf("add=");print_16(vResult16);printf("\n");
    printf("add=");print_32(vResult32);printf("\n");
    printf("add=");print_64(vResult64);printf("\n");

    vResult8 = _mm_sub_epi8(vOne8, vOne8);
    vResult16 = _mm_sub_epi16(vOne16, vOne16);
    vResult32 = _mm_sub_epi32(vOne32, vOne32);
    vResult64 = _mm_sub_epi64(vOne64, vOne64);

    printf("sub=");print_8(vResult8);printf("\n");
    printf("sub=");print_16(vResult16);printf("\n");
    printf("sub=");print_32(vResult32);printf("\n");
    printf("sub=");print_64(vResult64);printf("\n");

    vResult8 = _mm_cmpeq_epi8(vOne8, vOne8);
    vResult16 = _mm_cmpeq_epi16(vOne16, vOne16);
    vResult32 = _mm_cmpeq_epi32(vOne32, vOne32);
    vResult64 = _mm_cmpeq_epi64(vOne64, vOne64);

    printf("cmpeq 1 1 =");print_8(vResult8);printf("\n");
    printf("cmpeq 1 1 =");print_16(vResult16);printf("\n");
    printf("cmpeq 1 1 =");print_32(vResult32);printf("\n");
    printf("cmpeq 1 1 =");print_64(vResult64);printf("\n");

    vResult8 = _mm_cmpeq_epi8(vOne8, vZero8);
    vResult16 = _mm_cmpeq_epi16(vOne16, vZero16);
    vResult32 = _mm_cmpeq_epi32(vOne32, vZero32);
    vResult64 = _mm_cmpeq_epi64(vOne64, vZero64);

    printf("cmpeq 1 0 =");print_8(vResult8);printf("\n");
    printf("cmpeq 1 0 =");print_16(vResult16);printf("\n");
    printf("cmpeq 1 0 =");print_32(vResult32);printf("\n");
    printf("cmpeq 1 0 =");print_64(vResult64);printf("\n");

    vResult8 = _mm_cmpgt_epi8(vOne8, vOne8);
    vResult16 = _mm_cmpgt_epi16(vOne16, vOne16);
    vResult32 = _mm_cmpgt_epi32(vOne32, vOne32);
    vResult64 = _mm_cmpgt_epi64(vOne64, vOne64);

    printf("cmpgt 1 1 =");print_8(vResult8);printf("\n");
    printf("cmpgt 1 1 =");print_16(vResult16);printf("\n");
    printf("cmpgt 1 1 =");print_32(vResult32);printf("\n");
    printf("cmpgt 1 1 =");print_64(vResult64);printf("\n");

    vResult8 = _mm_cmpgt_epi8(vOne8, vZero8);
    vResult16 = _mm_cmpgt_epi16(vOne16, vZero16);
    vResult32 = _mm_cmpgt_epi32(vOne32, vZero32);
    vResult64 = _mm_cmpgt_epi64(vOne64, vZero64);

    printf("cmpgt 1 0 =");print_8(vResult8);printf("\n");
    printf("cmpgt 1 0 =");print_16(vResult16);printf("\n");
    printf("cmpgt 1 0 =");print_32(vResult32);printf("\n");
    printf("cmpgt 1 0 =");print_64(vResult64);printf("\n");

    vResult8 = _mm_cmplt_epi8(vOne8, vOne8);
    vResult16 = _mm_cmplt_epi16(vOne16, vOne16);
    vResult32 = _mm_cmplt_epi32(vOne32, vOne32);
    vResult64 = _mm_cmplt_epi64(vOne64, vOne64);

    printf("cmplt 1 1 =");print_8(vResult8);printf("\n");
    printf("cmplt 1 1 =");print_16(vResult16);printf("\n");
    printf("cmplt 1 1 =");print_32(vResult32);printf("\n");
    printf("cmplt 1 1 =");print_64(vResult64);printf("\n");

    vResult8 = _mm_cmplt_epi8(vOne8, vZero8);
    vResult16 = _mm_cmplt_epi16(vOne16, vZero16);
    vResult32 = _mm_cmplt_epi32(vOne32, vZero32);
    vResult64 = _mm_cmplt_epi64(vOne64, vZero64);

    printf("cmplt 1 0 =");print_8(vResult8);printf("\n");
    printf("cmplt 1 0 =");print_16(vResult16);printf("\n");
    printf("cmplt 1 0 =");print_32(vResult32);printf("\n");
    printf("cmplt 1 0 =");print_64(vResult64);printf("\n");

    printf("_mm_extract_epi8(vEnum,  0)=%d\n", _mm_extract_epi8(vEnum8,  0));
    printf("_mm_extract_epi8(vEnum,  1)=%d\n", _mm_extract_epi8(vEnum8,  1));
    printf("_mm_extract_epi8(vEnum,  2)=%d\n", _mm_extract_epi8(vEnum8,  2));
    printf("_mm_extract_epi8(vEnum,  3)=%d\n", _mm_extract_epi8(vEnum8,  3));
    printf("_mm_extract_epi8(vEnum,  4)=%d\n", _mm_extract_epi8(vEnum8,  4));
    printf("_mm_extract_epi8(vEnum,  5)=%d\n", _mm_extract_epi8(vEnum8,  5));
    printf("_mm_extract_epi8(vEnum,  6)=%d\n", _mm_extract_epi8(vEnum8,  6));
    printf("_mm_extract_epi8(vEnum,  7)=%d\n", _mm_extract_epi8(vEnum8,  7));
    printf("_mm_extract_epi8(vEnum,  8)=%d\n", _mm_extract_epi8(vEnum8,  8));
    printf("_mm_extract_epi8(vEnum,  9)=%d\n", _mm_extract_epi8(vEnum8,  9));
    printf("_mm_extract_epi8(vEnum, 10)=%d\n", _mm_extract_epi8(vEnum8, 10));
    printf("_mm_extract_epi8(vEnum, 11)=%d\n", _mm_extract_epi8(vEnum8, 11));
    printf("_mm_extract_epi8(vEnum, 12)=%d\n", _mm_extract_epi8(vEnum8, 12));
    printf("_mm_extract_epi8(vEnum, 13)=%d\n", _mm_extract_epi8(vEnum8, 13));
    printf("_mm_extract_epi8(vEnum, 14)=%d\n", _mm_extract_epi8(vEnum8, 14));
    printf("_mm_extract_epi8(vEnum, 15)=%d\n", _mm_extract_epi8(vEnum8, 15));

    printf("_mm_extract_epi16(vEnum, 0)=%d\n", _mm_extract_epi16(vEnum16, 0));
    printf("_mm_extract_epi16(vEnum, 1)=%d\n", _mm_extract_epi16(vEnum16, 1));
    printf("_mm_extract_epi16(vEnum, 2)=%d\n", _mm_extract_epi16(vEnum16, 2));
    printf("_mm_extract_epi16(vEnum, 3)=%d\n", _mm_extract_epi16(vEnum16, 3));
    printf("_mm_extract_epi16(vEnum, 4)=%d\n", _mm_extract_epi16(vEnum16, 4));
    printf("_mm_extract_epi16(vEnum, 5)=%d\n", _mm_extract_epi16(vEnum16, 5));
    printf("_mm_extract_epi16(vEnum, 6)=%d\n", _mm_extract_epi16(vEnum16, 6));
    printf("_mm_extract_epi16(vEnum, 7)=%d\n", _mm_extract_epi16(vEnum16, 7));

    printf("_mm_extract_epi32(vEnum, 0)=%d\n", _mm_extract_epi32(vEnum32, 0));
    printf("_mm_extract_epi32(vEnum, 1)=%d\n", _mm_extract_epi32(vEnum32, 1));
    printf("_mm_extract_epi32(vEnum, 2)=%d\n", _mm_extract_epi32(vEnum32, 2));
    printf("_mm_extract_epi32(vEnum, 3)=%d\n", _mm_extract_epi32(vEnum32, 3));

    printf("_mm_extract_epi64(vEnum, 0)=%lld\n", _mm_extract_epi64(vEnum64, 0));
    printf("_mm_extract_epi64(vEnum, 1)=%lld\n", _mm_extract_epi64(vEnum64, 1));

    printf("_mm_insert_epi8(vEnum, 99,  0)="); print_8(_mm_insert_epi8(vEnum8, 99,  0)); printf("\n");
    printf("_mm_insert_epi8(vEnum, 99,  1)="); print_8(_mm_insert_epi8(vEnum8, 99,  1)); printf("\n");
    printf("_mm_insert_epi8(vEnum, 99,  2)="); print_8(_mm_insert_epi8(vEnum8, 99,  2)); printf("\n");
    printf("_mm_insert_epi8(vEnum, 99,  3)="); print_8(_mm_insert_epi8(vEnum8, 99,  3)); printf("\n");
    printf("_mm_insert_epi8(vEnum, 99,  4)="); print_8(_mm_insert_epi8(vEnum8, 99,  4)); printf("\n");
    printf("_mm_insert_epi8(vEnum, 99,  5)="); print_8(_mm_insert_epi8(vEnum8, 99,  5)); printf("\n");
    printf("_mm_insert_epi8(vEnum, 99,  6)="); print_8(_mm_insert_epi8(vEnum8, 99,  6)); printf("\n");
    printf("_mm_insert_epi8(vEnum, 99,  7)="); print_8(_mm_insert_epi8(vEnum8, 99,  7)); printf("\n");
    printf("_mm_insert_epi8(vEnum, 99,  8)="); print_8(_mm_insert_epi8(vEnum8, 99,  8)); printf("\n");
    printf("_mm_insert_epi8(vEnum, 99,  9)="); print_8(_mm_insert_epi8(vEnum8, 99,  9)); printf("\n");
    printf("_mm_insert_epi8(vEnum, 99, 10)="); print_8(_mm_insert_epi8(vEnum8, 99, 10)); printf("\n");
    printf("_mm_insert_epi8(vEnum, 99, 11)="); print_8(_mm_insert_epi8(vEnum8, 99, 11)); printf("\n");
    printf("_mm_insert_epi8(vEnum, 99, 12)="); print_8(_mm_insert_epi8(vEnum8, 99, 12)); printf("\n");
    printf("_mm_insert_epi8(vEnum, 99, 13)="); print_8(_mm_insert_epi8(vEnum8, 99, 13)); printf("\n");
    printf("_mm_insert_epi8(vEnum, 99, 14)="); print_8(_mm_insert_epi8(vEnum8, 99, 14)); printf("\n");
    printf("_mm_insert_epi8(vEnum, 99, 15)="); print_8(_mm_insert_epi8(vEnum8, 99, 15)); printf("\n");

    printf("_mm_insert_epi16(vEnum, 99,  0)="); print_16(_mm_insert_epi16(vEnum16, 99,  0)); printf("\n");
    printf("_mm_insert_epi16(vEnum, 99,  1)="); print_16(_mm_insert_epi16(vEnum16, 99,  1)); printf("\n");
    printf("_mm_insert_epi16(vEnum, 99,  2)="); print_16(_mm_insert_epi16(vEnum16, 99,  2)); printf("\n");
    printf("_mm_insert_epi16(vEnum, 99,  3)="); print_16(_mm_insert_epi16(vEnum16, 99,  3)); printf("\n");
    printf("_mm_insert_epi16(vEnum, 99,  4)="); print_16(_mm_insert_epi16(vEnum16, 99,  4)); printf("\n");
    printf("_mm_insert_epi16(vEnum, 99,  5)="); print_16(_mm_insert_epi16(vEnum16, 99,  5)); printf("\n");
    printf("_mm_insert_epi16(vEnum, 99,  6)="); print_16(_mm_insert_epi16(vEnum16, 99,  6)); printf("\n");
    printf("_mm_insert_epi16(vEnum, 99,  7)="); print_16(_mm_insert_epi16(vEnum16, 99,  7)); printf("\n");

    printf("_mm_insert_epi32(vEnum, 99,  0)="); print_32(_mm_insert_epi32(vEnum32, 99,  0)); printf("\n");
    printf("_mm_insert_epi32(vEnum, 99,  1)="); print_32(_mm_insert_epi32(vEnum32, 99,  1)); printf("\n");
    printf("_mm_insert_epi32(vEnum, 99,  2)="); print_32(_mm_insert_epi32(vEnum32, 99,  2)); printf("\n");
    printf("_mm_insert_epi32(vEnum, 99,  3)="); print_32(_mm_insert_epi32(vEnum32, 99,  3)); printf("\n");

    printf("_mm_insert_epi64(vEnum, 99,  0)="); print_64(_mm_insert_epi64(vEnum64, 99,  0)); printf("\n");
    printf("_mm_insert_epi64(vEnum, 99,  1)="); print_64(_mm_insert_epi64(vEnum64, 99,  1)); printf("\n");

    _mm_store_si128(&vPointer[0], vEnum8);
    _mm_store_si128(&vPointer[1], vEnum8b);
    vResult8 = _mm_load_si128(&vPointer[0]);
    printf("load8a="); print_8(vResult8); printf("\n");
    vResult8 = _mm_load_si128(&vPointer[1]);
    printf("load8b="); print_8(vResult8); printf("\n");

    _mm_store_si128(&vPointer[0], vEnum16);
    _mm_store_si128(&vPointer[1], vEnum16b);
    vResult16 = _mm_load_si128(&vPointer[0]);
    printf("load16a="); print_16(vResult16); printf("\n");
    vResult16 = _mm_load_si128(&vPointer[1]);
    printf("load16b="); print_16(vResult16); printf("\n");

    _mm_store_si128(&vPointer[0], vEnum32);
    _mm_store_si128(&vPointer[1], vEnum32b);
    vResult32 = _mm_load_si128(&vPointer[0]);
    printf("load32a="); print_32(vResult32); printf("\n");
    vResult32 = _mm_load_si128(&vPointer[1]);
    printf("load32b="); print_32(vResult32); printf("\n");

    _mm_store_si128(&vPointer[0], vEnum64);
    _mm_store_si128(&vPointer[1], vEnum64b);
    vResult64 = _mm_load_si128(&vPointer[0]);
    printf("load64a="); print_64(vResult64); printf("\n");
    vResult64 = _mm_load_si128(&vPointer[1]);
    printf("load64b="); print_64(vResult64); printf("\n");

    printf("_mm_max_epi8(vOne, vZero)="); print_8(_mm_max_epi8(vOne8, vZero8)); printf("\n");
    printf("_mm_max_epi16(vOne, vZero)="); print_16(_mm_max_epi16(vOne16, vZero16)); printf("\n");
    printf("_mm_max_epi32(vOne, vZero)="); print_32(_mm_max_epi32(vOne32, vZero32)); printf("\n");
    printf("_mm_max_epi64(vOne, vZero)="); print_64(_mm_max_epi64(vOne64, vZero64)); printf("\n");

    printf("_mm_min_epi8(vOne, vZero)="); print_8(_mm_min_epi8(vOne8, vZero8)); printf("\n");
    printf("_mm_min_epi16(vOne, vZero)="); print_16(_mm_min_epi16(vOne16, vZero16)); printf("\n");
    printf("_mm_min_epi32(vOne, vZero)="); print_32(_mm_min_epi32(vOne32, vZero32)); printf("\n");
    printf("_mm_min_epi64(vOne, vZero)="); print_64(_mm_min_epi64(vOne64, vZero64)); printf("\n");


    vResult8 = _mm_cmpgt_epi8(vOne8, vOne8);
    vResult16 = _mm_cmpgt_epi16(vOne16, vOne16);
    vResult32 = _mm_cmpgt_epi32(vOne32, vOne32);
    vResult64 = _mm_cmpgt_epi64(vOne64, vOne64);

    printf("movemask cmpgt8 1 1 = %d\n", _mm_movemask_epi8(vResult8));
    printf("movemask cmpgt16 1 1 = %d\n", _mm_movemask_epi8(vResult16));
    printf("movemask cmpgt32 1 1 = %d\n", _mm_movemask_epi8(vResult32));
    printf("movemask cmpgt64 1 1 = %d\n", _mm_movemask_epi8(vResult64));

    vResult8 = _mm_cmpgt_epi8(vOne8, vZero8);
    vResult16 = _mm_cmpgt_epi16(vOne16, vZero16);
    vResult32 = _mm_cmpgt_epi32(vOne32, vZero32);
    vResult64 = _mm_cmpgt_epi64(vOne64, vZero64);

    printf("movemask cmpgt8 1 0 = %d\n", _mm_movemask_epi8(vResult8));
    printf("movemask cmpgt16 1 0 = %d\n", _mm_movemask_epi8(vResult16));
    printf("movemask cmpgt32 1 0 = %d\n", _mm_movemask_epi8(vResult32));
    printf("movemask cmpgt64 1 0 = %d\n", _mm_movemask_epi8(vResult64));

    vResult8 = _mm_and_si128(vResult8, vResult8);
    vResult16 = _mm_and_si128(vResult16, vResult16);
    vResult32 = _mm_and_si128(vResult32, vResult32);
    vResult64 = _mm_and_si128(vResult64, vResult64);

    printf("true and true ="); print_8(vResult8); printf("\n");
    printf("true and true ="); print_16(vResult16); printf("\n");
    printf("true and true ="); print_32(vResult32); printf("\n");
    printf("true and true ="); print_64(vResult64); printf("\n");

    vResult8 = _mm_or_si128(vResult8, vResult8);
    vResult16 = _mm_or_si128(vResult16, vResult16);
    vResult32 = _mm_or_si128(vResult32, vResult32);
    vResult64 = _mm_or_si128(vResult64, vResult64);

    printf("true or true ="); print_8(vResult8); printf("\n");
    printf("true or true ="); print_16(vResult16); printf("\n");
    printf("true or true ="); print_32(vResult32); printf("\n");
    printf("true or true ="); print_64(vResult64); printf("\n");

    vResult8 = _mm_andnot_si128(vResult8, vResult8);
    vResult16 = _mm_andnot_si128(vResult16, vResult16);
    vResult32 = _mm_andnot_si128(vResult32, vResult32);
    vResult64 = _mm_andnot_si128(vResult64, vResult64);

    printf("true andnot true ="); print_8(vResult8); printf("\n");
    printf("true andnot true ="); print_16(vResult16); printf("\n");
    printf("true andnot true ="); print_32(vResult32); printf("\n");
    printf("true andnot true ="); print_64(vResult64); printf("\n");

    vResult8 = _mm_xor_si128(vResult8, vResult8);
    vResult16 = _mm_xor_si128(vResult16, vResult16);
    vResult32 = _mm_xor_si128(vResult32, vResult32);
    vResult64 = _mm_xor_si128(vResult64, vResult64);

    printf("false xor false ="); print_8(vResult8); printf("\n");
    printf("false xor false ="); print_16(vResult16); printf("\n");
    printf("false xor false ="); print_32(vResult32); printf("\n");
    printf("false xor false ="); print_64(vResult64); printf("\n");

    vResult8 = _mm_slli_si128(vEnum8, 1);
    vResult16 = _mm_slli_si128(vEnum16, 2);
    vResult32 = _mm_slli_si128(vEnum32, 4);
    vResult64 = _mm_slli_si128(vEnum64, 8);

    printf("shift left enum ="); print_8(vResult8); printf("\n");
    printf("shift left enum ="); print_16(vResult16); printf("\n");
    printf("shift left enum ="); print_32(vResult32); printf("\n");
    printf("shift left enum ="); print_64(vResult64); printf("\n");

    vResult8 = _mm_rlli_si128(vEnum8, 1);
    vResult16 = _mm_rlli_si128(vEnum16, 2);
    vResult32 = _mm_rlli_si128(vEnum32, 4);
    vResult64 = _mm_rlli_si128(vEnum64, 8);

    printf("rotate left enum ="); print_8(vResult8); printf("\n");
    printf("rotate left enum ="); print_16(vResult16); printf("\n");
    printf("rotate left enum ="); print_32(vResult32); printf("\n");
    printf("rotate left enum ="); print_64(vResult64); printf("\n");

    vResult8 = _mm_srli_si128(vEnum8, 1);
    vResult16 = _mm_srli_si128(vEnum16, 2);
    vResult32 = _mm_srli_si128(vEnum32, 4);
    vResult64 = _mm_srli_si128(vEnum64, 8);

    printf("shift right enum ="); print_8(vResult8); printf("\n");
    printf("shift right enum ="); print_16(vResult16); printf("\n");
    printf("shift right enum ="); print_32(vResult32); printf("\n");
    printf("shift right enum ="); print_64(vResult64); printf("\n");

    printf("hmax enum = %d\n", _mm_hmax_epi8(vEnum8));
    printf("hmax enum = %d\n", _mm_hmax_epi16(vEnum16));
    printf("hmax enum = %d\n", _mm_hmax_epi32(vEnum32));
    printf("hmax enum = %lld\n", _mm_hmax_epi64(vEnum64));

    vMask8 = _mm_cmpeq_epi8(vZero8, vEnum8);
    vMask16 = _mm_cmpeq_epi16(vZero16, vEnum16);
    vMask32 = _mm_cmpeq_epi32(vZero32, vEnum32);
    vMask64 = _mm_cmpeq_epi64(vZero64, vEnum64);

    printf("mask cmpeq 0 enum ="); print_8(vMask8); printf("\n");
    printf("mask cmpeq 0 enum ="); print_16(vMask16); printf("\n");
    printf("mask cmpeq 0 enum ="); print_32(vMask32); printf("\n");
    printf("mask cmpeq 0 enum ="); print_64(vMask64); printf("\n");

    vResult8 = _mm_blendv_epi8(vOne8, vOne8, vMask8);
    vResult16 = _mm_blendv_epi8(vOne16, vOne16, vMask16);
    vResult32 = _mm_blendv_epi8(vOne32, vOne32, vMask32);
    vResult64 = _mm_blendv_epi8(vOne64, vOne64, vMask64);

    printf("blend vEnum vOne vMask ="); print_8(vResult8); printf("\n");
    printf("blend vEnum vOne vMask ="); print_16(vResult16); printf("\n");
    printf("blend vEnum vOne vMask ="); print_32(vResult32); printf("\n");
    printf("blend vEnum vOne vMask ="); print_64(vResult64); printf("\n");

    vResult8 = _mm_andnot_si128(vMask8, vOne8);
    vResult16 = _mm_andnot_si128(vMask16, vOne16);
    vResult32 = _mm_andnot_si128(vMask32, vOne32);
    vResult64 = _mm_andnot_si128(vMask64, vOne64);

    printf("andnot vMask vOne ="); print_8(vResult8); printf("\n");
    printf("andnot vMask vOne ="); print_16(vResult16); printf("\n");
    printf("andnot vMask vOne ="); print_32(vResult32); printf("\n");
    printf("andnot vMask vOne ="); print_64(vResult64); printf("\n");

    vResult8 = _mm_and_si128(vMask8, vOne8);
    vResult16 = _mm_and_si128(vMask16, vOne16);
    vResult32 = _mm_and_si128(vMask32, vOne32);
    vResult64 = _mm_and_si128(vMask64, vOne64);

    printf("and vMask vOne ="); print_8(vResult8); printf("\n");
    printf("and vMask vOne ="); print_16(vResult16); printf("\n");
    printf("and vMask vOne ="); print_32(vResult32); printf("\n");
    printf("and vMask vOne ="); print_64(vResult64); printf("\n");

#ifdef __POWER7__
    printf("__POWER7__ defined\n");
#endif
#ifdef __POWER8__
    printf("__POWER8__ defined\n");
#endif

    return 0;
}
