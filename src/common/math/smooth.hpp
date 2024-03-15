//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <cmath>

namespace math {

template<typename T>
static T MedianOf3(T u, T v, T w) {
    /* Median(u,v,w): */
    if ((u <= v && v <= w) ||
        (u >= v && v >= w))
        return v;
    if ((u <= w && w <= v) ||
        (u >= w && w >= v))
        return w;

    /* else */ return u;
}

/* Return (Index-1) of  median(u,v,w) , i.e.,
-1 : u
0 : v
1 : w
*/
template<typename T>
static int IndexOfMedianOf3(T u, T v, T w) {
    if ((u <= v && v <= w) ||
        (u >= v && v >= w))
        return 0;
    if ((u <= w && w <= v) ||
        (u >= w && w >= v))
        return 1;

    /* else */ return -1;
}

enum class SmoothEndRule {
    No,
    Copy,
    Tukey
};

template<typename T>
static bool SmoothEndStep(const T *x, T *y, size_t n, SmoothEndRule end_rule) {
    switch (end_rule) {
        default:
        case SmoothEndRule::No:
            return false;
        case SmoothEndRule::Copy:
            y[0] = x[0];
            y[n - 1] = x[n - 1];
            return false;
        case SmoothEndRule::Tukey: {
            bool chg = false;
            y[0] = MedianOf3(3 * y[1] - 2 * y[2], x[0], y[1]);
            chg = chg || (y[0] != x[0]);
            y[n - 1] = MedianOf3(y[n - 2], x[n - 1], 3 * y[n - 2] - 2 * y[n - 3]);
            chg = chg || (y[n - 1] != x[n - 1]);
            return chg;
        }
    }
}

template<typename T>
static bool Smooth3(const T *x, T *y, size_t n, SmoothEndRule end_rule) {
    // y[] := Running Median of three (x) = "3 (x[])" with "copy ends"
    // ---  return chg := ( y != x )
    bool chg = false;

    for (size_t i = 1; i < n - 1; i++) {
        int j = IndexOfMedianOf3(x[i - 1], x[i], x[i + 1]);
        y[i] = x[(int) i + j];
        chg = chg || j;
    }

    chg |= SmoothEndStep(x, y, n, end_rule);

    return chg;
}

template<typename T>
static size_t Smooth3R(const T *x, T *y, T *z, size_t n, SmoothEndRule end_rule) {
    // y[] := "3R"(x) ; 3R = Median of three, repeated until convergence
    size_t iter;
    bool chg;

    iter = chg = Smooth3(x, y, n, SmoothEndRule::Copy);

    while (chg) {
        if ((chg = Smooth3(y, z, n, SmoothEndRule::No))) {
            iter += 1;
            for (size_t i = 1; i < n - 1; i++)
                y[i] = z[i];
        }
    }

    chg |= SmoothEndStep(x, y, n, end_rule);

    return (iter ? iter : chg);
    /* = 0   <==>  only one "3" w/o any change
       = 1   <==>  either ["3" w/o change + endchange]
       or   [two "3"s, 2nd w/o change  ] */
}

template<typename T>
static bool SplitTest(const T *x, size_t i) {
    // Split test:
    //  Are we at a /-\ or \_/ location => split should be made ?

    if (x[i] != x[i + 1])
        return false;

    if ((x[i - 1] <= x[i] && x[i + 1] <= x[i + 2]) ||
        (x[i - 1] >= x[i] && x[i + 1] >= x[i + 2]))
        return false;

    /* else */ return true;
}

template<typename T>
static bool SmoothSplit3(const T *x, T *y, size_t n, bool do_ends) {
    // y[] := S(x[])  where S() = "sm_split3"
    bool chg = false;

    for (size_t i = 0; i < n; i++)
        y[i] = x[i];

    if (do_ends && SplitTest(x, 1)) {
        chg = true;
        y[1] = x[0];
        y[2] = MedianOf3(x[2], x[3], 3 * x[3] - 2 * x[4]);
    }

    for (size_t i = 2; i < n - 3; i++) {
        if (SplitTest(x, i)) {
            int j;
            // plateau at x[i] == x[i+1]

            // at left:
            if (-1 < (j = IndexOfMedianOf3(x[i], x[i - 1], 3 * x[i - 1] - 2 * x[i - 2]))) {
                y[i] = (j == 0 ? x[i - 1] : 3 * x[i - 1] - 2 * x[i - 2]);
                chg = (y[i] != x[i]);
            }

            // at right:
            if (-1 < (j = IndexOfMedianOf3(x[i + 1], x[i + 2], 3 * x[i + 2] - 2 * x[i + 3]))) {
                y[i + 1] = (j == 0 ? x[i + 2] : 3 * x[i + 2] - 2 * x[i + 3]);
                chg = (y[i + 1] != x[i + 1]);
            }
        }
    }

    if (do_ends && SplitTest(x, n - 3)) {
        chg = true;
        y[n - 2] = x[n - 1];
        y[n - 3] = MedianOf3(x[n - 3], x[n - 4], 3 * x[n - 4] - 2 * x[n - 5]);
    }

    return chg;
}

template<typename T>
size_t Smooth3RS3R(std::vector <T> &y, const std::vector <T> &x,
                   SmoothEndRule end_rule = SmoothEndRule::Tukey, bool split_ends = false) {
    // y[1:n] := "3R S 3R"(x[1:n]);  z = "work";
    size_t iter;
    bool chg;
    size_t n = x.size();

    y.resize(n);
    std::vector <T> z(n), w(n);

    iter = Smooth3R(&x[0], &y[0], &z[0], n, end_rule);
    chg = SmoothSplit3(&y[0], &z[0], n, split_ends);
    if (chg)
        iter += Smooth3R(&z[0], &y[0], &w[0], n, end_rule);

    /* else y == z already */
    return (iter + chg);
}

}
