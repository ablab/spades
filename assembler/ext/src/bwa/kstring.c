#include <stdarg.h>
#include <stdio.h>
#include "kstring.h"

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

#ifdef KSTRING_MAIN
#include <stdio.h>
int main()
{
	kstring_t *s;
	s = (kstring_t*)calloc(1, sizeof(kstring_t));
	ksprintf(s, "abcdefg: %d", 100);
	printf("%s\n", s->s);
	free(s);
	return 0;
}
#endif
