///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PARSE_SET_H
#define PARSE_SET_H

#include "CoreTools.h"

/// This file defines functions for parsing lists of entities, which are either
/// to be sorted or not.  There are three functions at present, ParseStringSet,
/// ParseIntSet, and ParseDoubleSet.
///
/// ParseStringSet: this is a MINIMAL implementation which could be enlarged.  At
/// present, it expects either:
/// - s, where s is any single string not starting with a curly bracket, or
/// - x{s1,...,sn}y [where x and y are possibly empty strings, not containing {}*],
///   which is expanded into the vector of strings xs1y,...,xsny
///   (in the same order), or
/// - {s1,...,sn}*k, same but expand to k copies of s1, k copies of s2, etc.
/// Also, any of s1,...,sn can have an individual multiplier, e.g., si*ki, which
/// expands to ki copies of just that element.
///
/// ParseIntSet: this defines a syntax for converting a string into a vec of 
/// integers, which is by default sorted, but can be left unsorted 
/// (by setting nosort = True).
///
/// In the current syntax there are several allowed formats:
///
/// [1] n, where n is an integer >= 0, represents a single integer
///
/// [2] {n1,...,nk}, where k >= 0 and n1,...,nk >= 0 are integers, representing the
///     given integers
///
/// [3] [a,b], where 0 <= a <= b, representing all integers between a and b, 
///     inclusive
///
/// [4] [a,b), where 0 <= a <= b, representing all integers x, a <= x < b
///
/// [5] @fn, were fn is a filename, representing the list of integers in file fn
///
/// [6] a combination of any of the above, separated by "|", representing the union
///     of the respective sets.
///
/// [7] random:n, where n is a nonnegative integer; this only makes sense if 
///     ParseIntSet is provided with start and stop arguments.
///
/// If start and stop are provided, values are checked to see that they lie on the
/// half open interval [start, stop).
///
/// This is not intended for applications where efficiency of parsing would matter.
///
/// NEGATIVE NUMBERS AND OVERFLOW ARE NOT CORRECTLY HANDLED AT PRESENT!
///
/// If status argument is provided, return even upon failure:
///      status = 0 returned upon success
///      status = 1 returned upon failure.
///
/// Note: if you use forms [2], [3], [4], or [6] to pass an argument to an 
/// executable on the command-line, you need to double-quote it.

void ParseIntSet( String descrip, vec<int>& answer, bool sortAnswer = true,
     const int start = 0, const int stop = -1 );

void ParseLongLongSet( String descrip, vec<longlong>& answer, 
     bool sortAnswer = true, const longlong start = 0, const longlong stop = -1 );

void ParseIntSet( String descrip, vec<int>& answer, int& status, 
     Bool ABORT_IF_BAD = False, bool sortAnswer = true, const int start = 0,
     const int stop = -1 );

void ParseLongLongSet( String descrip, vec<longlong>& answer, int& status, 
     Bool ABORT_IF_BAD = False, bool sortAnswer = true, const longlong start = 0,
     const longlong stop = -1 );

void ParseStringSet( String descrip, vec<String>& answer, bool recurse = false );

/// Read doubles from a String and put them in a vec.
/// The allowed syntaxes (as described above) are
/// [1] x, where x is a string convertible to a double
/// [2] {x1,x2,...,xk}, where k >= 0 and x1,...,xk are strings convertible to doubles
/// and
/// [5] @fn, where fn is the name of a file containing whitespace-separated strings 
///     convertible to doubles

void ParseDoubleSet( 
     const String & descrip, vec<double> & answer, bool sortAnswer = true );

#endif
