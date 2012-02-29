///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file Bases.cc
 * \author tsharpe
 * \date Apr 10, 2009
 *
 * \brief Classes representing nucleotides.
 */
#include "dna/Bases.h"
#include "random/Random.h"
#include "system/System.h"

GeneralizedBase const GeneralizedBase::R('R','Y',true,5,"AG");
GeneralizedBase const GeneralizedBase::Y('Y','R',true,10,"CT");
GeneralizedBase const GeneralizedBase::K('K','M',true,12,"GT");
GeneralizedBase const GeneralizedBase::M('M','K',true,3,"AC");
GeneralizedBase const GeneralizedBase::S('S','S',true,6,"CG");
GeneralizedBase const GeneralizedBase::W('W','W',true,9,"AT");
GeneralizedBase const GeneralizedBase::B('B','V',true,14,"CGT");
GeneralizedBase const GeneralizedBase::V('V','B',true,7,"ACG");
GeneralizedBase const GeneralizedBase::D('D','H',true,13,"AGT");
GeneralizedBase const GeneralizedBase::H('H','D',true,11,"ACT");
GeneralizedBase const GeneralizedBase::N('N','N',true,15,"ACGT");
GeneralizedBase const GeneralizedBase::X('-','-',true,0,"");

Base const Base::A('A','T',false,1,"A",BASE_A);
Base const Base::C('C','G',false,2,"C",BASE_C);
Base const Base::G('G','C',false,4,"G",BASE_G);
Base const Base::T('T','A',false,8,"T",BASE_T);

#define BA Base::A
#define BC Base::C
#define BG Base::G
#define BT Base::T

GeneralizedBase const*const GeneralizedBase::gCharToGenBase[256] =
{
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //00
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //10
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &X, &N,  0, //20
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //30
     0,&BA, &B,&BC, &D,  0,  0,&BG, &H,  0,  0, &K,  0, &M, &N,  0, //40
     0,  0, &R, &S,&BT,  0, &V, &W, &X, &Y,  0,  0,  0,  0,  0,  0, //50
     0,&BA, &B,&BC, &D,  0,  0,&BG, &H,  0,  0, &K,  0, &M, &N,  0, //60
     0,  0, &R, &S,&BT,&BT, &V, &W, &X, &Y,  0,  0,  0,  0,  0,  0, //70
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //80
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //90
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //a0
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //b0
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //c0
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //d0
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //e0
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0  //f0
};

GeneralizedBase const*const GeneralizedBase::gBitsToGenBase[16] =
{ &X, &BA, &BC, &M, &BG, &R, &S, &V, &BT, &W, &Y, &H, &K, &D, &B, &N };

char GeneralizedBase::gRevComp[256] =
{
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //00
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //10
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,'.','-',  0, //20
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //30
     0,'T','V','G','H',  0,  0,'C','D',  0,  0,'M',  0,'K','N',  0, //40
     0,  0,'Y','S','A',  0,'B','W','X','R',  0,  0,  0,  0,  0,  0, //50
     0,'t','v','g','h',  0,  0,'c','d',  0,  0,'m',  0,'k','n',  0, //60
     0,  0,'y','s','a','t','b','w','x','r',  0,  0,  0,  0,  0,  0, //70
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //80
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //90
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //a0
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //b0
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //c0
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //d0
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //e0
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0  //f0
};

Base const*const Base::gValToBase[4] =
{ &A, &C, &G, &T };

Base const*const Base::gCharToBase[256] =
{
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //00
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //10
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //20
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //30
     0, &A,  0, &C,  0,  0,  0, &G,  0,  0,  0,  0,  0,  0,  0,  0, //40
     0,  0,  0,  0, &T, &T,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //50
     0, &A,  0, &C,  0,  0,  0, &G,  0,  0,  0,  0,  0,  0,  0,  0, //60
     0,  0,  0,  0, &T, &T,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //70
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //80
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //90
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //a0
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //b0
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //c0
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //d0
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //e0
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0  //f0
};

unsigned char const Base::gRCByte[256] = {
    0xff, 0xbf, 0x7f, 0x3f, 0xef, 0xaf, 0x6f, 0x2f,
    0xdf, 0x9f, 0x5f, 0x1f, 0xcf, 0x8f, 0x4f, 0x0f,
    0xfb, 0xbb, 0x7b, 0x3b, 0xeb, 0xab, 0x6b, 0x2b,
    0xdb, 0x9b, 0x5b, 0x1b, 0xcb, 0x8b, 0x4b, 0x0b,
    0xf7, 0xb7, 0x77, 0x37, 0xe7, 0xa7, 0x67, 0x27,
    0xd7, 0x97, 0x57, 0x17, 0xc7, 0x87, 0x47, 0x07,
    0xf3, 0xb3, 0x73, 0x33, 0xe3, 0xa3, 0x63, 0x23,
    0xd3, 0x93, 0x53, 0x13, 0xc3, 0x83, 0x43, 0x03,
    0xfe, 0xbe, 0x7e, 0x3e, 0xee, 0xae, 0x6e, 0x2e,
    0xde, 0x9e, 0x5e, 0x1e, 0xce, 0x8e, 0x4e, 0x0e,
    0xfa, 0xba, 0x7a, 0x3a, 0xea, 0xaa, 0x6a, 0x2a,
    0xda, 0x9a, 0x5a, 0x1a, 0xca, 0x8a, 0x4a, 0x0a,
    0xf6, 0xb6, 0x76, 0x36, 0xe6, 0xa6, 0x66, 0x26,
    0xd6, 0x96, 0x56, 0x16, 0xc6, 0x86, 0x46, 0x06,
    0xf2, 0xb2, 0x72, 0x32, 0xe2, 0xa2, 0x62, 0x22,
    0xd2, 0x92, 0x52, 0x12, 0xc2, 0x82, 0x42, 0x02,
    0xfd, 0xbd, 0x7d, 0x3d, 0xed, 0xad, 0x6d, 0x2d,
    0xdd, 0x9d, 0x5d, 0x1d, 0xcd, 0x8d, 0x4d, 0x0d,
    0xf9, 0xb9, 0x79, 0x39, 0xe9, 0xa9, 0x69, 0x29,
    0xd9, 0x99, 0x59, 0x19, 0xc9, 0x89, 0x49, 0x09,
    0xf5, 0xb5, 0x75, 0x35, 0xe5, 0xa5, 0x65, 0x25,
    0xd5, 0x95, 0x55, 0x15, 0xc5, 0x85, 0x45, 0x05,
    0xf1, 0xb1, 0x71, 0x31, 0xe1, 0xa1, 0x61, 0x21,
    0xd1, 0x91, 0x51, 0x11, 0xc1, 0x81, 0x41, 0x01,
    0xfc, 0xbc, 0x7c, 0x3c, 0xec, 0xac, 0x6c, 0x2c,
    0xdc, 0x9c, 0x5c, 0x1c, 0xcc, 0x8c, 0x4c, 0x0c,
    0xf8, 0xb8, 0x78, 0x38, 0xe8, 0xa8, 0x68, 0x28,
    0xd8, 0x98, 0x58, 0x18, 0xc8, 0x88, 0x48, 0x08,
    0xf4, 0xb4, 0x74, 0x34, 0xe4, 0xa4, 0x64, 0x24,
    0xd4, 0x94, 0x54, 0x14, 0xc4, 0x84, 0x44, 0x04,
    0xf0, 0xb0, 0x70, 0x30, 0xe0, 0xa0, 0x60, 0x20,
    0xd0, 0x90, 0x50, 0x10, 0xc0, 0x80, 0x40, 0x00
};

unsigned char GeneralizedBase::chooseRandom() const
{
    if ( this == &X )
        FatalErr("Can't select a random base from a 'no base' call.");
    return Base::char2Val(mBases[randomx() % (mEndBases - mBases)]);
}

void GeneralizedBase::fromBitsFatalErr( unsigned char bits )
{
    FatalErr("Incorrect bit value for generalized base: " << static_cast<unsigned int>(bits));
}

void GeneralizedBase::fromCharFatalErr( char chr )
{
    FatalErr("Incorrect character value for generalized base: " << chr
        << " [int value: " << int(chr) << "]");
}

void Base::fromValFatalErr( unsigned char val )
{
    FatalErr("Incorrect code for base: " << static_cast<unsigned int>(val)
         << " = \'" << val << "\'" );
}

void Base::fromCharFatalErr( char chr )
{
    FatalErr("Incorrect character value for base: " << chr);
}
