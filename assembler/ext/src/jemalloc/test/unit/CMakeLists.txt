###########################################################################
# Build Unit tests
# This links to exported C_JETLIB static version of the library

set(UnitTests
    SFMT.c
    a0.c
    arena_reset.c
    arena_reset_prof.c
    atomic.c
    background_thread.c
    background_thread_enable.c
    base.c
    binshard.c
    bit_util.c
    bitmap.c
    ckh.c
    decay.c
    div.c
    emitter.c
    extent_quantize.c
    extent_util.c
    fork.c
    hash.c
    hook.c
    huge.c
    junk.c
    junk_alloc.c
    junk_free.c
    log.c
    mallctl.c
    malloc_io.c
    math.c
    mq.c
    mtx.c
    nstime.c
    pack.c
    pages.c
    ph.c
    prng.c
    prof_accum.c
    prof_active.c
    prof_gdump.c
    prof_idump.c
    prof_log.c
    prof_reset.c
    prof_tctx.c
    prof_thread_name.c
    ql.c
    qr.c
    rb.c
    retained.c
    rtree.c
    sc.c
    seq.c
    size_classes.c
    slab.c
    smoothstep.c
    spin.c
    stats.c
    stats_print.c
    test_hooks.c
    ticker.c
    tsd.c
    witness.c
    zero.c
    )

foreach(sourcefile ${UnitTests})
    get_filename_component(exename ${sourcefile} NAME_WE)
    set(exename "${exename}${install_suffix}")
    createTest(${exename} ${sourcefile} JEMALLOC_UNIT_TEST)
    target_link_libraries(${PROJECT_NAME}-test-${exename} ${PROJECT_NAME}UnitTest)
endforeach()

target_link_libraries(${PROJECT_NAME}-test-math m)
