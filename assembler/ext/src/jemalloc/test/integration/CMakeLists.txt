set(TESTS_INTEGRATION
    MALLOCX_ARENA.c
    aligned_alloc.c
    allocated.c
    extent.c
    malloc.c
    mallocx.c
    overflow.c
    posix_memalign.c
    rallocx.c
    sdallocx.c
    slab_sizes.c
    thread_arena.c
    thread_tcache_enabled.c
    xallocx.c
    )

if(JEMALLOC_EXPERIMENTAL_SMALLOCX_API)
    list(APPEND TESTS_INTEGRATION smallocx.c)
endif()

##########################################################################
# Build integration tests
# This links to SO library and UTIL_INTEGRATION
foreach(sourcefile ${TESTS_INTEGRATION})
    get_filename_component(exename ${sourcefile} NAME_WE)
    set(exename "${exename}${install_suffix}")
    createTest(${exename} ${sourcefile} JEMALLOC_INTEGRATION_TEST)
endforeach(sourcefile ${TESTS_UNIT})
