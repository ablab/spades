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
    bool use_name     : 1; // false implies use_quality = false, use_comment = false, true implies use_comment = true
    bool use_comment  : 1; // can override use_comment = false in case when use_name = false
    bool use_quality  : 1;
    bool validate     : 1;
    bool paired       : 1;

    static FileReadFlags empty() {
        return { PhredOffset,
                 /* use name */ false, /* use comment */ false,
                 /* use quality */ false, /* validate */ false };
    }

    static FileReadFlags only_names() {
        return { PhredOffset,
                 /* use name */ true, /* use comment */ false,
                 /* use quality */ false, /* validate */ false };
    }

    static FileReadFlags only_comments() {
        return { PhredOffset,
                 /* use name */ false, /* use comment */ true,
                 /* use quality */ false, /* validate */ false };
    }

    static FileReadFlags names_and_comments() {
        return { PhredOffset,
                 /* use name */ true, /* use comment */ true,
                 /* use quality */ false, /* validate */ false };
    }

    
    FileReadFlags()
            : offset(PhredOffset),
              use_name(true), use_comment(true), use_quality(true), validate(true), paired(false) {}
    FileReadFlags(OffsetType o)
            : offset(o),
              use_name(true), use_comment(true), use_quality(true), paired(false) {}
    FileReadFlags(OffsetType o, bool n, bool q)
            : offset(o), use_name(n), use_comment(n), use_quality(q), paired(false) {}
    FileReadFlags(OffsetType o, bool n, bool c, bool q, bool v)
            : offset(o), use_name(n), use_comment(c), use_quality(q), validate(v), paired(false) {}
};

}
