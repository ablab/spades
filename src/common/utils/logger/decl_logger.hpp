//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2020-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

inline const char* __scope_source_name() {
    return " General ";
}

#define DECL_LOGGER(source)                                         \
static const char* __scope_source_name() {                          \
    return source;                                                  \
}
