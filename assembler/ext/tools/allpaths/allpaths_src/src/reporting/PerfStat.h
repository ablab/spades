///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file PerfStat.h
 * \author tsharpe
 * \date Jul 21, 2010
 *
 * \brief A measurement of performance to log.
 */
#ifndef REPORTING_PERFSTAT_H_
#define REPORTING_PERFSTAT_H_

#include <iostream>
#include <string>

/// A measurement of performance.
///
/// The only real utility is to record these to a log, which will be digested
/// by some external tool.  In particular, Dexter can run commands under control
/// of a performance tester that will read the log and ascertain that the
/// performance is up to snuff.  We do this before generating nightly packages
/// to try to prevent unleashing broken code upon the world.
///
/// Use it like this:
/// PerfStat::log() << PerfStat("FrobDeg","degree of frobnication",frobDeg());
///
/// The name can be any token without whitespace (this is not checked).  You
/// might wish to scope these in some sensible way so that sub-programs and
/// sub-sub-programs don't overwrite each others' results:  i.e., the name
/// "AllPathsLG.Contig.N50" is probably a better name than "N50".
/// The gloss is just a human readable description of the measurement.  It
/// shouldn't contain newline characters (not checked).
///
/// If your value is very large or small, you may want to use manipulators
/// to set fixed or scientific format, precision, etc.
/// E.g.,
/// PerfStat::log() << std::fixed << std::setprecision(2)
///     << PerfStat("PI","current value of pi",3.14);
///
/// You can read a performance statistics log using an ifstream and operator>>,
/// but I don't know why you'd want to.
///
class PerfStat
{
    typedef std::istream istream;
    typedef std::ostream ostream;
    typedef std::string string;

public:
    PerfStat() : mVal(0.) {}

    PerfStat( string const& name, string const& gloss, long double val )
    : mName(name), mGloss(gloss), mVal(val) {}

    // compiler-supplied copying and destructor is OK

    friend istream& operator>>( istream& is, PerfStat& perfStat );
    friend ostream& operator<<( ostream& os, PerfStat const& perfStat );

    static ostream& log();

private:
    string mName;
    string mGloss;
    long double mVal;
};



// usage:
//
// PerfStat::log() << PerfStatBlockStart("some random table");
// ...
// cout << "whatever" << endl;
// ...
// 
// PerfStat::log() << PerfStatBlockStop();
//
std::string PerfStatBlockStart(const std::string & block_name);
std::string PerfStatBlockStop();









#endif /* REPORTING_PERFSTAT_H_ */
