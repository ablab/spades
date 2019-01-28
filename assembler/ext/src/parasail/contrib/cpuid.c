/* https://github.com/Mysticial/Flops/blob/e571da6e94f7b6d2d1a90e87b19398c5c4de4375/version1/source/cpuid.c */

/* cpuid.c
 * 
 * Author           : Alexander J. Yee
 * Date Created     : 01/19/2012
 * Last Modified    : 01/25/2012
 * 
 * 
 * 
 * And of course... The typical copyright stuff...
 * 
 *      Redistribution of this program in both source or binary, regardless of
 *  form, with or without modification is permitted as long as the following
 *  conditions are met:
 *          1.  This copyright notice is maintained either inline in the source
 *              or distributed with the binary.
 *          2.  A list of all contributing authors along with their contributions
 *              is included either inline in the source or distributed with the
 *              binary.
 *          3.  The following disclaimer is maintained either inline in the
 *              source or distributed with the binary.
 * 
 *  Disclaimer:
 *  This software is provided "as is", without any guarantee made to its
 *  suitability or fitness for any particular use. It may contain bugs so use
 *  of this program is at your own risk. I take no responsibility for any
 *  damage that may unintentionally be caused through its use.
 */

#ifndef _cpuid_c
#define _cpuid_c
#include <stdio.h>
#include "cpuid.h"
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
#ifdef WIN32
#else
void cpuid(int *info,int x){
    int ax,bx,cx,dx;

    __asm__ __volatile__ ("cpuid": "=a" (ax), "=b" (bx), "=c" (cx), "=d" (dx) : "a" (x));

    info[0] = ax;
    info[1] = bx;
    info[2] = cx;
    info[3] = dx;
}
#endif
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void cpuid_print_name(){
    int name[13];
    cpuid(name + 0,0x80000002);
    cpuid(name + 4,0x80000003);
    cpuid(name + 8,0x80000004);
    name[12] = '\0';

    printf("CPU Name = %s\n",(char*)name);
    printf("\n");
}
void cpuid_print_exts(){
    int x64     = 0;
    int MMX     = 0;
    int SSE     = 0;
    int SSE2    = 0;
    int SSE3    = 0;
    int SSSE3   = 0;
    int SSE41   = 0;
    int SSE42   = 0;
    int SSE4a   = 0;
    int AVX     = 0;
    int XOP     = 0;
    int FMA3    = 0;
    int FMA4    = 0;

    int info[4];
    cpuid(info, 0);
    int nIds = info[0];

    cpuid(info, 0x80000000);
    int nExIds = info[0];

    //  Detect Instruction Set
    if (nIds >= 1){
        cpuid(info,0x00000001);
        MMX   = (info[3] & ((int)1 << 23)) != 0;
        SSE   = (info[3] & ((int)1 << 25)) != 0;
        SSE2  = (info[3] & ((int)1 << 26)) != 0;
        SSE3  = (info[2] & ((int)1 <<  0)) != 0;

        SSSE3 = (info[2] & ((int)1 <<  9)) != 0;
        SSE41 = (info[2] & ((int)1 << 19)) != 0;
        SSE42 = (info[2] & ((int)1 << 20)) != 0;

        AVX   = (info[2] & ((int)1 << 28)) != 0;
        FMA3  = (info[2] & ((int)1 << 12)) != 0;
    }

    if (nExIds >= 0x80000001){
        cpuid(info,0x80000001);
        x64   = (info[3] & ((int)1 << 29)) != 0;
        SSE4a = (info[2] & ((int)1 <<  6)) != 0;
        FMA4  = (info[2] & ((int)1 << 16)) != 0;
        XOP   = (info[2] & ((int)1 << 11)) != 0;
    }

    printf("Hardware Features:\n");
    printf("x64   = %d\n",x64);
    printf("MMX   = %d\n",MMX);
    printf("SSE   = %d\n",SSE);
    printf("SSE2  = %d\n",SSE2);
    printf("SSE3  = %d\n",SSE3);
    printf("SSSE3 = %d\n",SSSE3);
    printf("SSE4a = %d\n",SSE4a);
    printf("SSE41 = %d\n",SSE41);
    printf("SSE42 = %d\n",SSE42);
    printf("AVX   = %d\n",AVX);
    printf("FMA3  = %d\n",FMA3);
    printf("FMA4  = %d\n",FMA4);
    printf("XOP   = %d\n",XOP);
    printf("\n");
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
#endif

