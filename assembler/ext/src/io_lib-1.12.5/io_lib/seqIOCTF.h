#ifndef _seqIOCTF_h_
#define _seqIOCTF_h_

#include "io_lib/Read.h"

#ifdef __cplusplus
extern "C" {
#endif

Read *ctfFRead (mFILE *ff) ;
int ctfFWrite (mFILE *ff, Read *read) ;

#ifdef __cplusplus
}
#endif

#endif /* _seqIOCTF_h_ */















