/**
 * @file
 *
 * This file exists solely to prevent errors with
 * macOS's ar when producing empty archives:
 *
 *   rm -f src/libnovec_plain_sse2.a && ar csr src/libnovec_plain_sse2.a
 *   ar: no archive members specified
 *
 * Copyright (c) 2017 Battelle Memorial Institute.
 */

static void dummy_function_for_macos_ar(void) {}
