############################################################################
# Common source for Unit, Integration and stress test libraries
set(UnitTestCommon
    SFMT.c
    btalloc.c
    btalloc_0.c
    btalloc_1.c
    math.c
    mq.c
    mtx.c
    test.c
    thd.c
    timer.c
    )

############################################################################
# Unit tests test library
# Unit test have access to all of the interfaces
add_library(${PROJECT_NAME}UnitTest STATIC ${UnitTestCommon})
add_dependencies(${PROJECT_NAME}UnitTest ${PROJECT_NAME}-jet-static)
target_include_directories(${PROJECT_NAME}UnitTest PRIVATE
    ${PROJECT_SOURCE_DIR}/test/include)
target_compile_definitions(${PROJECT_NAME}UnitTest
    PRIVATE
    JEMALLOC_UNIT_TEST)

############################################################################
# Integration test library
add_library(${PROJECT_NAME}IntegrationTest STATIC ${UnitTestCommon})
add_dependencies(${PROJECT_NAME}IntegrationTest ${PROJECT_NAME}-jet-static)
target_include_directories(${PROJECT_NAME}IntegrationTest PRIVATE
    ${PROJECT_SOURCE_DIR}/test/include)
target_compile_definitions(${PROJECT_NAME}IntegrationTest PRIVATE
    JEMALLOC_INTEGRATION_TEST)

############################################################################
# Stress test library to link with microbench
add_library(${PROJECT_NAME}TestStress STATIC ${UnitTestCommon})
add_dependencies(${PROJECT_NAME}TestStress ${PROJECT_NAME}-jet-static)

target_include_directories(${PROJECT_NAME}TestStress PRIVATE
    ${PROJECT_SOURCE_DIR}/test/include)
target_compile_definitions(${PROJECT_NAME}TestStress PRIVATE
    JEMALLOC_STRESS_TEST
    JEMALLOC_STRESS_TESTLIB)
