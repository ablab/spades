//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * nucmer_coords_tool.hpp
 *
 *  Created on: 22.08.2011

 */

#ifndef NUCMER_COORDS_TOOL_HPP_
#define NUCMER_COORDS_TOOL_HPP_
#include <stdio.h>
#include <cstdlib>
#include <string>
#include <vector>
#include <math.h>

struct NucmerAllign {
    int inReferenceStart;
    int inReferenceEnd;
    int inContigStart;
    int inContigEnd;
    int inReferenceLength;
    int inContigLength;
    double ReferenceQuality;
    std::string contigName;
};
std::vector<NucmerAllign> ReadCoordsFile(std::string file_name);
inline bool MyLess(NucmerAllign i, NucmerAllign j){
    return (i.inReferenceStart < j.inReferenceStart);
}

inline bool MyLess2(NucmerAllign i, NucmerAllign j){
    if (i.contigName < j.contigName) return true;
    if (i.contigName > j.contigName) return false;
    if (i.inContigLength > j.inContigLength) return true;
    if (i.inContigLength < j.inContigLength) return false;
    if (i.ReferenceQuality > j.ReferenceQuality) return true;
    if (i.ReferenceQuality < j.ReferenceQuality) return false;
    return (i.inReferenceStart < j.inReferenceStart);
}
inline bool Redundant(NucmerAllign i, NucmerAllign j){
    if (i.contigName != j.contigName) return false;
    int i_min = std::min(i.inContigStart, i.inContigEnd);
    int i_max = std::max(i.inContigStart, i.inContigEnd);
    int j_min = std::min(j.inContigStart, j.inContigEnd);
    int j_max = std::max(j.inContigStart, j.inContigEnd);
    if ((i_min <= j_min)&&(i_max >= j_max)&&(i.ReferenceQuality >= j.ReferenceQuality)) return true;
    return false;

}


#endif /* NUCMER_COORDS_TOOL_HPP_ */
