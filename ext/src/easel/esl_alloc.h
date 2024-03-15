#ifndef eslALLOC_INCLUDED
#define eslALLOC_INCLUDED

#include <stdlib.h>

extern void *esl_alloc_aligned(size_t size, size_t alignment);
extern void  esl_alloc_free(void *p);

#endif // eslALLOC_INCLUDED
