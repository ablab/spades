##############################################################################
## Determine system page size and cache line size at compile time
function (GetSystemAddrBits OUTPUT_VAR_VADDR)
    # Direct all the files into one folder
    set(WORK_FOLDER "${PROJECT_BINARY_DIR}/GetSystemAddrBits")
    file(MAKE_DIRECTORY ${WORK_FOLDER})

    set(SRC "${WORK_FOLDER}/getsystemaddrbits.c")
    set(COMPILE_OUTPUT_FILE "${WORK_FOLDER}/getsystemaddrbits.log")

    file(WRITE ${SRC} "
    #include <stdio.h>
#ifdef _WIN32
    #include <limits.h>
    #include <intrin.h>
    typedef unsigned __int32 uint32_t;
#else
    #include <stdint.h>
#endif

int main(void) {
    uint32_t r[4];
    uint32_t eax_in = 0x80000008U;
#ifdef _WIN32
    __cpuid((int *)r, (int)eax_in);
#else
    asm volatile (\"cpuid\"
        : \"=a\" (r[0]), \"=b\" (r[1]), \"=c\" (r[2]), \"=d\" (r[3])
        : \"a\" (eax_in), \"c\" (0)
    );
#endif

    uint32_t eax_out = r[0];
    uint32_t vaddr = ((eax_out & 0x0000ff00U) >> 8);
    printf(\"%u\", vaddr);
    return 0;
}
")

    try_run(RUN_RESULT COMPILE_RESULT
        "${WORK_FOLDER}"
        "${SRC}"
        COMPILE_OUTPUT_VARIABLE COMPILE_OUTPUT
        RUN_OUTPUT_VARIABLE RUN_OUTPUT
        )

    if(NOT COMPILE_RESULT)
        file(WRITE ${COMPILE_OUTPUT_FILE} ${COMPILE_OUTPUT})
        message(FATAL_ERROR "GetSystemAddrBits failed compilation see ${COMPILE_OUTPUT_FILE}")
    endif()

    if("${RUN_RESULT}" STREQUAL "FAILED_TO_RUN")
        message(FATAL_ERROR "GetSystemAddrBits failed to run executable")
    endif()

    message(STATUS "System addr bits: ${RUN_OUTPUT}")

    set(${OUTPUT_VAR_VADDR} ${RUN_OUTPUT} PARENT_SCOPE)
endfunction(GetSystemAddrBits)
