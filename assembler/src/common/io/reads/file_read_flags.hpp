//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

namespace io {

/*
* This enumerate contains offset type.
* UnknownOffset is equal to "offset = 0".
* PhredOffset is equal to "offset = 33".
* SolexaOffset is equal to "offset = 64".
*/
//todo change to enum class
enum OffsetType {
    UnknownOffset = 0,
    PhredOffset = 33,
    SolexaOffset = 64
};

struct FileReadFlags {
    OffsetType offset : 8;
    bool use_name     : 1;
    bool use_quality  : 1;
    bool validate     : 1;

    FileReadFlags()
            : offset(PhredOffset), use_name(true), use_quality(true), validate(true) {}
    FileReadFlags(OffsetType o)
            : offset(o), use_name(true), use_quality(true) {}
    FileReadFlags(OffsetType o, bool n, bool q)
            : offset(o), use_name(n), use_quality(q) {}
    FileReadFlags(OffsetType o, bool n, bool q, bool v)
            : offset(o), use_name(n), use_quality(q), validate(v) {}

};

}
