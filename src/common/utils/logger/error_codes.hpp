//***************************************************************************
//* Copyright (c) 2025 SPAdes team
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <cstdint>

// Predefined error codes range from 64 to 127 will not trigger message suggesting to report a bug
// Predefined error codes range from 128 to 255 (system calls, negative values) and
// from 1 to 63 will trigger message suggesting to report a bug
enum ErrorCodes: std::int8_t {
    InvalidInputFormat = 64,
    InputFileNotFound = 65,
    IOError = 66,
    InvalidParameter = 67,
    MemoryLimitExceeded = 68
};

