///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file PerfStat.cc
 * \author tsharpe
 * \date Jul 21, 2010
 *
 * \brief
 */
#include "reporting/PerfStat.h"
#include <fstream>
#include <stdlib.h>
#include <string.h>

using std::ios_base;


std::istream& operator>>( std::istream& is, PerfStat& perfStat )
{
    is >> perfStat.mName;
    is >> perfStat.mVal;
    is >> std::ws;
    getline(is,perfStat.mGloss);
    return is;
}

std::ostream& operator<<( std::ostream& os, PerfStat const& perfStat )
{
  os << "PERFSTAT: " << perfStat.mGloss << " [" << perfStat.mName << "] = " << perfStat.mVal << std::endl;
  return os;
}

std::ostream& PerfStat::log()
{
  static int dfltPrecision = std::cout.precision();
  std::cout.precision(dfltPrecision);
  std::cout.unsetf(ios_base::floatfield);
  return std::cout;
}



std::string PerfStatBlockStart(const std::string & block_name)
{
  return "PERFSTAT: BLOCK_START [" + block_name + "]\n";
}

std::string PerfStatBlockStop()
{
  return "PERFSTAT: BLOCK_STOP\n";
}

