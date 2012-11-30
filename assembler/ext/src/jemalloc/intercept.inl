#include <jemalloc/jemalloc.h>

#ifdef HAVE_FEATURES_H
# include <features.h>
#endif
#ifdef HAVE_SYS_CDEFS_H
# include <sys/cdefs.h.>
#endif

#if defined(__APPLE__)
/* Nothing to do here, everything will be initialized via ctor attribute inside zone.c */
#elif defined(__GNUC__)
# define ALIAS(fn) __attribute__((alias(#fn)))
void* malloc(size_t size) __THROW               ALIAS(je_malloc);
void free(void* ptr) __THROW                    ALIAS(je_free);
void* realloc(void* ptr, size_t size) __THROW   ALIAS(je_realloc);
void* calloc(size_t n, size_t size) __THROW     ALIAS(je_calloc);
void cfree(void* ptr) __THROW                   ALIAS(je_free);
void* memalign(size_t align, size_t s) __THROW  ALIAS(je_memalign);
void* valloc(size_t size) __THROW               ALIAS(je_valloc);
void* pvalloc(size_t size) __THROW              ALIAS(je_valloc);
int posix_memalign(void** r, size_t a, size_t s) __THROW ALIAS(je_posix_memalign);

# if defined(__GLIBC__)
void* __libc_malloc(size_t size)                      ALIAS(je_malloc);
void  __libc_free(void* ptr)                          ALIAS(je_free);
void* __libc_realloc(void* ptr, size_t size)          ALIAS(je_realloc);
void* __libc_calloc(size_t n, size_t size)            ALIAS(je_calloc);
void  __libc_cfree(void* ptr)                         ALIAS(je_free);
void* __libc_memalign(size_t align, size_t s)         ALIAS(je_memalign);
void* __libc_valloc(size_t size)                      ALIAS(je_valloc);
void* __libc_pvalloc(size_t size)                     ALIAS(je_valloc);
int   __posix_memalign(void** r, size_t a, size_t s)  ALIAS(je_posix_memalign);

#  include <malloc.h>
static void* glibc_override_malloc(size_t size, const void *caller) {
  return je_malloc(size);
}
static void* glibc_override_realloc(void *ptr, size_t size,
                                           const void *caller) {
  return je_realloc(ptr, size);
}
static void glibc_override_free(void *ptr, const void *caller) {
  je_free(ptr);
}
static void* glibc_override_memalign(size_t align, size_t size,
                                     const void *caller) {
  return je_memalign(align, size);
}

/* From GNU libc 2.14 this macro is defined, to declare
   hook variables as volatile. Define it as empty for
   older glibc versions */
#ifndef __MALLOC_HOOK_VOLATILE
# define __MALLOC_HOOK_VOLATILE
#endif

void *(*__MALLOC_HOOK_VOLATILE __malloc_hook)(size_t size, const void *caller) = &glibc_override_malloc;
void *(*__MALLOC_HOOK_VOLATILE __realloc_hook)(void *ptr, size_t size, const void *caller) = &glibc_override_realloc;
void (*__MALLOC_HOOK_VOLATILE __free_hook) (void*, const void *) = &glibc_override_free;
void *(*__MALLOC_HOOK_VOLATILE __memalign_hook) (size_t alignment, size_t size, const void *caller) = &glibc_override_memalign;

# endif /* __GLIBC__ */
#undef ALIAS
#else
# error "Intercepting for this platform is not supported"
#endif
