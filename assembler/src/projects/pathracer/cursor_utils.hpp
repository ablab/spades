//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <parallel_hashmap/phmap.h>
#include <vector>

#include "common/utils/verify.hpp"

template <typename Cursor>
auto vertex_cursors(const std::vector<Cursor> &cursors, typename Cursor::Context context) {
    phmap::flat_hash_set<Cursor> result;

    for (const Cursor &cursor : cursors) {
        VERIFY(!cursor.is_empty());
        if (cursor.prev(context).size() != 1 || cursor.next(context).size() != 1) {
            result.insert(cursor);
        }
    }

    INFO("Isolated loops detection");
    // detect and break isolated loops
    phmap::flat_hash_set<Cursor> processed;
    for (const Cursor &cursor : cursors) {
        if (result.count(cursor) || processed.count(cursor)) {
            continue;
        }

        Cursor p = cursor;
        while (!result.count(p) && !processed.count(p)) {
            processed.insert(p);
            p = p.prev(context)[0];
        }
        if (p == cursor) {
            INFO("Isolated loop detected");
            result.insert(cursor);
        }
    }
    return result;
}

template <typename Cursor>
auto extract_leftmost_cursors(const phmap::flat_hash_set<Cursor> &cursors, const phmap::flat_hash_set<Cursor> &vertices,
                              typename Cursor::Context context) {
    // Exclude all trivial cursors having preceding cursor on the same edge in the input set
    std::vector<Cursor> result;

    for (const Cursor &cursor : cursors) {
        VERIFY(!cursor.is_empty());
        if (vertices.count(cursor)) {
            result.push_back(cursor);
            continue;
        }

        Cursor p = cursor;
        do {
            VERIFY(p.prev(context).size() == 1 && p.next(context).size() == 1);
            p = p.prev(context)[0];
        } while (!cursors.count(p) && !vertices.count(p));

        if (!cursors.count(p)) {
            result.push_back(cursor);
        }
    }

    return result;
}
