#.rst:
# Findxgetbv
# --------
#
# Finds xgetbv support
#
# This module can be used to detect xgetbv support in a C compiler.
#
# The following variables are set:
#
# ::
#
#    HAVE_XGETBV - true if xgetbv is detected
#
#=============================================================================

# sample xgetbv source code to test
set(XGETBV_C_TEST_SOURCE
"
#include <stdint.h>
int check_xcr0_ymm()
{
    uint32_t xcr0;
#if defined(_MSC_VER)
    xcr0 = (uint32_t)_xgetbv(0);  /* min VS2010 SP1 compiler is required */
#else
    __asm__ (\"xgetbv\" : \"=a\" (xcr0) : \"c\" (0) : \"%edx\" );
#endif
    return ((xcr0 & 6) == 6); /* checking if xmm and ymm state are enabled in XCR0 */
}
int main(void) { return check_xcr0_ymm(); }
")

check_c_source_compiles("${XGETBV_C_TEST_SOURCE}" HAVE_XGETBV)

