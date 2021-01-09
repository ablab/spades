#if !defined(MI_IN_ALLOC_C)
#error "this file should be included from 'alloc.c' (so aliases can work)"
#endif

#if defined(MI_MALLOC_OVERRIDE_WRAP)

#if defined(_WIN32) || defined(__MACH__)
#error "this file is supported only on Linux"
#endif

#if !defined(MI_FORWARD)
  #if (defined(__GNUC__) && __GNUC__ >= 9)
    #define MI_FORWARD(fun)      __attribute__((alias(#fun), used, visibility("default"), copy(fun)))
  #else
    #define MI_FORWARD(fun)      __attribute__((alias(#fun), used, visibility("default")))
  #endif
#endif

#if (defined(__GNUC__) || defined(__clang__)) && !defined(__MACH__)
#pragma GCC visibility push(default)
#endif

#ifdef __cplusplus
extern "C" {
#endif

void* __wrap_malloc(size_t size)              MI_FORWARD(mi_malloc);
void* __wrap_calloc(size_t size, size_t n)    MI_FORWARD(mi_calloc);
void* __wrap_realloc(void* p, size_t newsize) MI_FORWARD(mi_realloc);
void  __wrap_free(void* p)                    MI_FORWARD(mi_free);

void   __wrap_cfree(void* p)                    MI_FORWARD(mi_free);
void*  __wrap_reallocf(void* p, size_t newsize) MI_FORWARD(mi_reallocf);
size_t __wrap_malloc_size(const void* p)        MI_FORWARD(mi_usable_size);
size_t __wrap__malloc_usable_size(void *p)      MI_FORWARD(mi_usable_size);

void* __wrap_valloc(size_t size)                                     { return mi_valloc(size); }
void* __wrap_pvalloc(size_t size)                                    { return mi_pvalloc(size); }
void* __wrap_reallocarray(void* p, size_t count, size_t size)        { return mi_reallocarray(p, count, size); }
void* __wrap_memalign(size_t alignment, size_t size)                 { return mi_memalign(alignment, size); }
int   __wrap_posix_memalign(void** p, size_t alignment, size_t size) { return mi_posix_memalign(p, alignment, size); }
void* __wrap__aligned_malloc(size_t alignment, size_t size)          { return mi_aligned_alloc(alignment, size); }
void* __wrap_aligned_alloc(size_t alignment, size_t size)            { return mi_aligned_alloc(alignment, size); }

#if defined(__GLIBC__) && defined(__linux__)
#ifndef __MALLOC_HOOK_VOLATILE
#define __MALLOC_HOOK_VOLATILE /**/
#endif

// Setup malloc hooks. These should never be called, but still...
static void* glibc_override_malloc(size_t size, const void*) {
  return mi_malloc(size);
}
static void* glibc_override_realloc(void* ptr, size_t size,
                                    const void*) {
  return mi_realloc(ptr, size);
}
static void glibc_override_free(void* ptr, const void*) {
  mi_free(ptr);
}
static void* glibc_override_memalign(size_t align, size_t size,
                                     const void*) {
  return mi_memalign(align, size);
}

void* (*__MALLOC_HOOK_VOLATILE __malloc_hook)(size_t, const void*) =
    &glibc_override_malloc;
void* (*__MALLOC_HOOK_VOLATILE __realloc_hook)(void*, size_t, const void*) =
    &glibc_override_realloc;
void (*__MALLOC_HOOK_VOLATILE __free_hook)(void*,
                                           const void*) = &glibc_override_free;
void* (*__MALLOC_HOOK_VOLATILE __memalign_hook)(size_t, size_t, const void*) =
    &glibc_override_memalign;
#endif

#ifdef __cplusplus
}
#endif

#if (defined(__GNUC__) || defined(__clang__)) && !defined(__MACH__)
#pragma GCC visibility pop
#endif

#endif // defined(MI_MALLOC_OVERRIDE_WRAP)
