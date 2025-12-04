//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

// Predefined error codes range from 64 to 127 will not trigger message suggesting to report a bug
// Predefined error codes range from 128 to 255 (negative values) will trigger message suggesting to report a bug
enum ErrorCodes {
    InvalidInputFormat = 64,
    InputFileNotFound = 65,
    IOError = 66,
    InvalidParameter = 67,
    MemoryLimitExceeded = 68
};

