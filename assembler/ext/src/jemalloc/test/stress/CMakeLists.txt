############################################################################
# Microbench
set(TESTS_STRESS
    hookbench.c
    microbench.c
    )

# Build benchmark
# This consumes C_JETLIB, Link to SO
foreach(sourcefile ${TESTS_STRESS})
    get_filename_component(exename ${sourcefile} NAME_WE)
    set(exename "${exename}${install_suffix}")
    createTest(${exename} ${sourcefile} JEMALLOC_STRESS_TEST)
    target_link_libraries(${PROJECT_NAME}-test-${exename} ${PROJECT_NAME}TestStress)
endforeach()
